function traj = solve_scenario_bvp(scn)
% SOLVE_SCENARIO_BVP  Solve one test maneuver with the PMP / bvp4c method of optimizeBVP.m.
%
%   traj = solve_scenario_bvp(scn)
%
% Called by run_scenarios.m, which builds scn from its knobs. You shouldn't need to edit
% this file: the machine limits live in machine_limits() below and stay fixed.
%
% scn fields (poses are [trolley_m, cable_m, slew_deg], positions are ground-frame X/Y in m):
%   name          short name, used for the output files
%   start, finish rest poses
%   obstacles     struct array, each either
%                   kind 'wall': a, b (wall end points [X Y]), thickness, height
%                   kind 'pin' : xy ([X Y]), radius, height
%   waypoints     2 x K payload X/Y points the starting guess passes through ([] = none)
%   guess_shape   'arc' (trolley, cable, slew blended directly) or 'straight' (payload in a
%                 straight line); with waypoints the guess is always straight segments
%   tf_guess      starting guess for the move time in s ([] = estimate from the limits)
%   time_weight   speed vs. smoothness (3 = near minimum time, as in optimizeBVP.m)
%   preview_only  true: only draw the guess and the obstacles, don't solve
%   obstacle_stiffness_start  barrier weight for the first solves (1e-2 normally; 1 for a
%                 loop that starts and ends at the same pose, so it can't cut through)
%   out_dir       folder for <name>.mat and <name>.xlsx
%
% traj is the same format optimizeBVP.m saves (so sim_trajectory.m can play it), plus
% scenario info and metrics. It is saved to out_dir/<name>.mat and also copied to
% trajectory.mat next to this file. Empty if preview_only.

M = machine_limits();
here = fileparts(mfilename('fullpath'));

xInitial = pose_to_state(scn.start,  'start',  scn.name, M);
xFinal   = pose_to_state(scn.finish, 'finish', scn.name, M);

[obs, obs_xy, obs_rho] = build_obstacles(scn.obstacles, M);
check_endpoint_clearance(xInitial, xFinal, obs, M, scn.name);

%% Starting guess
N_init = 60;
tau = linspace(0, 1, N_init);
x_guess = make_guess(xInitial, xFinal, scn.waypoints, scn.guess_shape, tau);
if isempty(scn.tf_guess)
    tf_guess = estimate_tf(x_guess, M);
else
    tf_guess = scn.tf_guess;
end
x_guess(6:8, :) = x_guess(6:8, :) / tf_guess;          % make_guess returned d/dtau
guess_clear = check_guess_clearance(x_guess, obs_xy, obs_rho, scn.name);

if abs(x_guess(3, end) - xFinal(3)) > 1e-6
    warning(['%s: the guess path goes around the mast to slew %.0f deg, but the finish pose ' ...
             'says %.0f deg. Add or subtract 360 from the finish slew if that is not what you meant.'], ...
            scn.name, rad2deg(x_guess(3, end)), rad2deg(xFinal(3)));
end

if scn.preview_only
    preview_guess(scn, x_guess, obs, obs_xy, obs_rho, M);
    fprintf('%s: preview only. Guess t_f = %.1f s, guess clearance %.0f mm. Not solved.\n', ...
            scn.name, tf_guess, 1000*guess_clear);
    traj = [];
    return
end

%% Symbolic Hamiltonian (rebuilt only when the obstacles change)
[ode_rhs_fn, u_star_fn, H_val_fn] = pmp_functions(M, obs_xy, obs_rho);

bvp_opts = bvpset('RelTol', 1e-3, 'AbsTol', 1e-5, 'NMax', 15000, 'Stats', 'off');
ws = scn.time_weight / 3.0;               % optimizeBVP.m's schedule ends at w_t = 3

%% Stage 0: fixed t_f, finds the costates
fprintf('\n=== %s ===\n', scn.name);
fprintf('Stage 0 (fixed t_f = %.1f s, finding costates)...\n', tf_guess);
k_obs0 = scn.obstacle_stiffness_start;       % 1e-2 normally; stiffer keeps a loop from collapsing
p_stage0 = [0.20*ws, 1e-6, k_obs0];
Y_guess = [x_guess; zeros(10, N_init)];
sol0 = bvpinit(tau, @(t) interp1(tau, Y_guess', t)');
lastwarn('');
sol0 = bvp4c(@(t, Y) bvp_ode(t, Y, tf_guess, ode_rhs_fn, p_stage0), ...
             @(Ya, Yb) [Ya(1:10) - xInitial; Yb(1:10) - xFinal], sol0, bvp_opts);
if ~isempty(lastwarn)
    error(['%s: the first solve did not converge (%s).\nTry: a different tf_guess, a larger ' ...
           'guess_clearance, or moving the start/finish poses away from the obstacles.'], scn.name, lastwarn);
end
%% Stage 0b: free t_f
% Freeing t_f straight from a poor guess sometimes fails. Then: for a fixed-t_f optimum,
% dJ/dt_f = H(t_f), so H > 0 means a shorter move is cheaper and H < 0 a longer one. Step
% the fixed t_f (warm-starting each solve from the last) until H changes sign, and free
% t_f from the closer of the two.
free_tf = @(s0, tf0) try_free_tf(s0, tf0, ode_rhs_fn, H_val_fn, xInitial, xFinal, p_stage0, bvp_opts);
sol = free_tf(sol0, tf_guess);
if isempty(sol)
    H_end = @(s) H_val_fn(s.y(1:10, end), s.y(11:20, end), p_stage0(1), p_stage0(2), p_stage0(3));
    best = sol0;  tf_best = tf_guess;  H_best = H_end(sol0);
    fixed_bc = @(Ya, Yb) [Ya(1:10) - xInitial; Yb(1:10) - xFinal];
    while true
        tf_try = tf_best * (1 + 0.15*sign(-H_best));
        if tf_try < 0.3*tf_guess || tf_try > 3*tf_guess, break; end
        lastwarn('');
        try
            s_try = bvp4c(@(t, Y) bvp_ode(t, Y, tf_try, ode_rhs_fn, p_stage0), fixed_bc, best, bvp_opts);
            ok = isempty(lastwarn);
        catch
            ok = false;
        end
        if ~ok, break; end
        H_try = H_end(s_try);
        crossed = sign(H_try) ~= sign(H_best);
        if ~crossed || abs(H_try) < abs(H_best)
            best = s_try;  tf_best = tf_try;  H_best = H_try;
        end
        if crossed, break; end
    end
    fprintf('Stage 0b: walked the fixed t_f to %.2f s (H(t_f) = %+.2e), freeing it from there...\n', tf_best, H_best);
    sol = free_tf(best, tf_best);
    if isempty(sol)
        error(['%s: could not free the move time.\nTry a different tf_guess or a larger ' ...
               'guess_clearance.'], scn.name);
    end
end

%% Continuation: [w_t, k_wall, k_obs], each stage reached in adaptive sub-steps
stages = [0.20*ws, 1e-5, 1e-2;
          0.50*ws, 1e-4, 1e-1;
          1.50*ws, 5e-4, 1.0;
          3.00*ws, 1e-3, 1.0;
          3.00*ws, 1.0,  1.0];
stages(:, 3) = max(stages(:, 3), k_obs0);
p_prev = p_stage0;
for s = 1:size(stages, 1)
    p_to = stages(s, :);
    a = 0;  step = 1;
    while a < 1
        a_try = min(1, a + step);
        p_try = blend(p_prev, p_to, a_try);
        lastwarn('');
        try
            sol_try = bvp4c(@(t, Y, tf) bvp_ode(t, Y, tf, ode_rhs_fn, p_try), ...
                @(Ya, Yb, tf) bvp_bc(Ya, Yb, xInitial, xFinal, H_val_fn, p_try), sol, bvp_opts);
            ok = isempty(lastwarn);
        catch err
            ok = false;  fprintf('    (%s)\n', err.message);
        end
        if ok
            sol = sol_try;  a = a_try;  step = min(2*step, 1);
        else
            step = step/2;
            fprintf('    step failed, halving to %.3f\n', step);
            if step < 1/64
                error(['%s: stuck in stage %d at %.0f%%.\nTry: a different tf_guess, a smaller ' ...
                       'time_weight, or a guess path farther from the obstacles.'], scn.name, s, 100*a);
            end
        end
    end
    p_prev = p_to;
    fprintf('Stage %d/%d done: t_f = %.2f s, %d mesh points\n', s, size(stages, 1), sol.parameters, numel(sol.x));
end

%% Unpack, check, save
tf_opt = sol.parameters;
t_grid = sol.x * tf_opt;
x_opt  = sol.y(1:10, :);
lam_opt = sol.y(11:20, :);
u_opt = zeros(3, numel(t_grid));
for k = 1:numel(t_grid), u_opt(:, k) = u_star_fn(x_opt(:, k), lam_opt(:, k)); end
H_t = arrayfun(@(k) H_val_fn(x_opt(:, k), lam_opt(:, k), p_prev(1), p_prev(2), p_prev(3)), 1:numel(t_grid));

Q = payload_xyz(x_opt, M.H_jib);
swing_deg = rad2deg(acos(min(1, cos(x_opt(4, :)).*cos(x_opt(5, :)))));
min_clear = inf;
for j = 1:numel(obs)
    d = arrayfun(@(i) obstacle_sdf(Q(:, i), obs(j), M.R_pay), 1:size(Q, 2));
    min_clear = min(min_clear, min(d));
end

fprintf('Converged: t_f = %.2f s   (PMP check: H(t) spread %.1e, should be ~0)\n', tf_opt, max(abs(H_t)));
fprintf('  peak speed / limit: trolley %.0f%%, hoist %.0f%%, slew %.0f%%\n', ...
        100*max(abs(x_opt(6, :)))/M.dx_max, 100*max(abs(x_opt(7, :)))/M.dl_max, 100*max(abs(x_opt(8, :)))/M.dth_max);
fprintf('  peak swing during the move: %.2f deg', max(swing_deg));
if isfinite(min_clear), fprintf(',  closest approach to an obstacle: %.0f mm', 1000*min_clear); end
fprintf('\n');

traj.t       = t_grid;
traj.x       = x_opt;
traj.u       = u_opt;
traj.lam     = lam_opt;
traj.tf      = tf_opt;
traj.w_t     = stages(end, 1);
traj.N       = numel(t_grid) - 1;
traj.x_lb    = [M.x_min; M.l_min; -inf; -inf; -inf; -M.dx_max; -M.dl_max; -M.dth_max; -inf; -inf];
traj.x_ub    = [M.x_max; M.l_max;  inf;  inf;  inf;  M.dx_max;  M.dl_max;  M.dth_max;  inf;  inf];
traj.u_lb    = -M.u_max;  traj.u_ub = M.u_max;
traj.obs     = obs;
traj.H_jib   = M.H_jib;
traj.R_pay   = M.R_pay;
traj.created = char(datetime('now'));
traj.scenario        = scn;
traj.peak_swing_deg  = max(swing_deg);
traj.min_clearance_m = min_clear;

if ~isfolder(scn.out_dir), mkdir(scn.out_dir); end
save(fullfile(scn.out_dir, [scn.name '.mat']), '-struct', 'traj');
save(fullfile(here, 'trajectory.mat'), '-struct', 'traj');          % what sim_trajectory.m plays
export_xlsx(fullfile(scn.out_dir, [scn.name '.xlsx']), t_grid, u_opt, M.u_max);
fprintf('  saved %s.mat / .xlsx in %s\n', scn.name, scn.out_dir);
end

%% ========================================================================
%  Machine (fixed)
%  ========================================================================
function M = machine_limits()
M.x_min = 0.0;  M.x_max = 0.925;                 % trolley range (m)
M.l_min = 0.0;  M.l_max = 1.450;                 % cable range (m)
M.dx_max  = 72.6e-3;                             % m/s
M.dl_max  = 30.0e-3;                             % m/s
M.dth_max = deg2rad(32.4);                       % rad/s
M.u_max   = [72.6e-3; 30.0e-3; deg2rad(32.4)];   % [ddx; ddl; ddth] max accelerations
M.H_jib = 2.0;                                   % jib height above the floor (m)
M.R_pay = 0.08;                                  % payload radius, used as clearance (m)
M.circle_spacing = 0.05;                         % walls become circles <= 5 cm apart
end

function x = pose_to_state(pose, which, name, M)
if numel(pose) ~= 3
    error('%s: the %s pose must be [trolley_m, cable_m, slew_deg].', name, which);
end
if pose(1) < M.x_min + 0.02 || pose(1) > M.x_max - 0.02
    error('%s: %s trolley = %.3f m is outside the usable range %.2f to %.3f m.', ...
          name, which, pose(1), M.x_min + 0.02, M.x_max - 0.02);
end
if pose(2) < 0.2 || pose(2) > M.l_max - 0.02
    error('%s: %s cable = %.3f m is outside the usable range 0.2 to %.2f m.', name, which, pose(2), M.l_max - 0.02);
end
x = [pose(1); pose(2); deg2rad(pose(3)); zeros(7, 1)];
end

%% ========================================================================
%  Obstacles
%  ========================================================================
function [obs, obs_xy, obs_rho] = build_obstacles(spec, M)
% obs: the real shapes (boxes/cylinders, for the sim and the clearance check).
% obs_xy / obs_rho: the smooth-barrier circles the PMP solver sees. They act as infinitely
% tall, so the solver keeps the payload out of an obstacle's footprint at every height.
obs = struct('type', {}, 'center', {}, 'half', {}, 'yaw', {}, 'xy', {}, 'R', {}, 'z', {});
obs_xy = zeros(2, 0);  obs_rho = zeros(1, 0);
for j = 1:numel(spec)
    o = spec(j);
    switch o.kind
        case 'wall'
            a = o.a(:);  b = o.b(:);  L = norm(b - a);
            if L < o.thickness, error('Wall %d: its two end points are closer than its thickness.', j); end
            d = (b - a) / L;
            obs(end+1).type = 'box'; %#ok<AGROW>
            obs(end).center = [(a + b)/2; o.height/2];
            obs(end).half   = [L/2; o.thickness/2; o.height/2];
            obs(end).yaw    = atan2(d(2), d(1));
            rho = o.thickness/2 + M.R_pay;
            e1 = a + d*o.thickness/2;  e2 = b - d*o.thickness/2;     % circles end flush with the wall ends
            n = max(2, ceil(norm(e2 - e1) / M.circle_spacing) + 1);
            obs_xy  = [obs_xy, e1 + (e2 - e1)*linspace(0, 1, n)]; %#ok<AGROW>
            obs_rho = [obs_rho, rho*ones(1, n)]; %#ok<AGROW>
        case 'pin'
            obs(end+1).type = 'cylinder'; %#ok<AGROW>
            obs(end).xy = o.xy(:);  obs(end).R = o.radius;  obs(end).z = [0; o.height];
            obs_xy  = [obs_xy, o.xy(:)]; %#ok<AGROW>
            obs_rho = [obs_rho, o.radius + M.R_pay]; %#ok<AGROW>
        otherwise
            error('Unknown obstacle kind ''%s''.', o.kind);
    end
end
end

function check_endpoint_clearance(xI, xF, obs, M, name)
names = {'start', 'finish'};  X = [xI, xF];
for k = 1:2
    p = payload_xyz(X(:, k), M.H_jib);
    for j = 1:numel(obs)
        d = obstacle_sdf(p, obs(j), M.R_pay);
        if d <= 0.005
            error('%s: the %s pose is inside (or touching) obstacle %d (clearance %.0f mm). Move the pose or the obstacle.', ...
                  name, names{k}, j, 1000*d);
        end
    end
end
end

function c = check_guess_clearance(xg, obs_xy, obs_rho, name)
c = inf;
if isempty(obs_rho), return; end
P = [xg(1, :).*cos(xg(3, :)); xg(1, :).*sin(xg(3, :))];
for j = 1:size(obs_xy, 2)
    c = min(c, min(vecnorm(P - obs_xy(:, j))) - obs_rho(j));
end
if c <= 0
    error(['%s: the starting guess passes through an obstacle (%.0f mm inside). Increase ' ...
           'guess_clearance, pick the other side, or move the poses. Set PREVIEW_ONLY = true to see it.'], ...
          name, -1000*c);
end
end

%% ========================================================================
%  Starting guess
%  ========================================================================
function xg = make_guess(xI, xF, wp, shape, tau)
% States along the guess path with a smoothstep time profile (at rest at both ends).
% Rates in rows 6:8 are d/dtau here; the caller divides by t_f.
s = 3*tau.^2 - 2*tau.^3;
xg = zeros(10, numel(tau));
if isempty(wp) && strcmp(shape, 'arc')
    for i = 1:3, xg(i, :) = xI(i) + (xF(i) - xI(i))*s; end
else
    P0 = xI(1)*[cos(xI(3)); sin(xI(3))];
    P1 = xF(1)*[cos(xF(3)); sin(xF(3))];
    P  = [P0, wp, P1];
    cum = [0, cumsum(vecnorm(diff(P, 1, 2)))];  cum = cum / cum(end);
    XY = interp1(cum, P.', s).';
    xg(1, :) = vecnorm(XY);
    th = unwrap(atan2(XY(2, :), XY(1, :)));
    xg(3, :) = th - th(1) + xI(3);
    xg(2, :) = xI(2) + (xF(2) - xI(2))*s;
    if min(xg(1, :)) < 0.03
        error('The straight-line guess passes over the mast. Use guess_shape = ''arc'' or add a waypoint.');
    end
end
for i = 1:3, xg(i+5, :) = gradient(xg(i, :), tau); end
end

function tf = estimate_tf(xg, M)
% Time for each drive to cover the guess path's travel at its limits (trapezoid profile),
% plus margin for settling the swing.
v = [M.dx_max; M.dl_max; M.dth_max];
t_axis = 0;
for i = 1:3
    d = sum(abs(diff(xg(i, :))));
    if d >= v(i)^2 / M.u_max(i), t = d/v(i) + v(i)/M.u_max(i);
    else,                        t = 2*sqrt(d / M.u_max(i));
    end
    t_axis = max(t_axis, t);
end
tf = 1.3*t_axis + 2.0;
end

function preview_guess(scn, xg, obs, obs_xy, obs_rho, M)
figure('Name', ['Preview: ' scn.name], 'Color', 'w');  hold on;  axis equal;  grid on;
a = linspace(0, 2*pi, 200);
plot(M.x_max*cos(a), M.x_max*sin(a), ':', 'Color', [.5 .5 .5]);
plot(0, 0, 'k^', 'MarkerFaceColor', 'k');
for j = 1:numel(obs)
    switch obs(j).type
        case 'box'
            [cx, cy] = box_corners(obs(j));  fill(cx, cy, [.85 .2 .2], 'FaceAlpha', .4, 'EdgeColor', [.85 .2 .2]);
        case 'cylinder'
            fill(obs(j).xy(1) + obs(j).R*cos(a), obs(j).xy(2) + obs(j).R*sin(a), [.85 .2 .2], 'FaceAlpha', .4);
    end
end
for j = 1:size(obs_xy, 2)
    plot(obs_xy(1, j) + obs_rho(j)*cos(a), obs_xy(2, j) + obs_rho(j)*sin(a), '-', 'Color', [1 .6 .6]);
end
P = [xg(1, :).*cos(xg(3, :)); xg(1, :).*sin(xg(3, :))];
plot(P(1, :), P(2, :), 'b-', 'LineWidth', 1.5);
if ~isempty(scn.waypoints), plot(scn.waypoints(1, :), scn.waypoints(2, :), 'bx', 'MarkerSize', 10); end
plot(P(1, 1), P(2, 1), 'go', 'MarkerFaceColor', 'g');
plot(P(1, end), P(2, end), 'ks', 'MarkerFaceColor', 'k');
xlabel('X (m)');  ylabel('Y (m)');
title({['Preview: ' strrep(scn.name, '_', ' ')], ...
       'blue = starting guess, pink = keep-out zone (payload center), dotted = trolley reach'});
end

function [cx, cy] = box_corners(ob)
c = [-1 1 1 -1; -1 -1 1 1] .* ob.half(1:2);
R = [cos(ob.yaw) -sin(ob.yaw); sin(ob.yaw) cos(ob.yaw)];
v = R*c + ob.center(1:2);
cx = v(1, :);  cy = v(2, :);
end

function Q = payload_xyz(X, H)
% Payload center in the ground frame, Z up from the floor.
x = X(1, :);  l = X(2, :);  th = X(3, :);  phi = X(4, :);  psi = X(5, :);
rx = x + l.*sin(phi).*cos(psi);  ry = l.*sin(psi);
Q = [rx.*cos(th) - ry.*sin(th);  rx.*sin(th) + ry.*cos(th);  H - l.*cos(phi).*cos(psi)];
end

%% ========================================================================
%  BVP pieces (as in optimizeBVP.m)
%  ========================================================================
function dY = bvp_ode(~, Y, tf, ode_rhs_fn, p)
dY = tf * ode_rhs_fn(Y(1:10), Y(11:20), p(1), p(2), p(3));
end

function res = bvp_bc(Ya, Yb, xI, xF, H_val_fn, p)
% 10 initial states + 10 final states + H(t_f) = 0 (free final time)
res = [Ya(1:10) - xI; Yb(1:10) - xF; H_val_fn(Yb(1:10), Yb(11:20), p(1), p(2), p(3))];
end

function sol = try_free_tf(s0, tf0, ode_rhs_fn, H_val_fn, xI, xF, p, opts)
% One free-t_f solve from a fixed-t_f solution; [] if bvp4c fails or misses its tolerance.
lastwarn('');
try
    sol = bvp4c(@(t, Y, tf) bvp_ode(t, Y, tf, ode_rhs_fn, p), ...
                @(Ya, Yb, tf) bvp_bc(Ya, Yb, xI, xF, H_val_fn, p), ...
                bvpinit(s0.x, @(t) deval(s0, t), tf0), opts);
    if ~isempty(lastwarn), sol = []; end
catch
    sol = [];
end
end

function p = blend(p_from, p_to, a)
% w_t linearly, the stiffness weights geometrically
p = p_from + a*(p_to - p_from);
for i = 2:3
    if p_from(i) > 0 && p_to(i) > 0, p(i) = p_from(i) * (p_to(i)/p_from(i))^a; end
end
end

function [ode_rhs_fn, u_star_fn, H_val_fn] = pmp_functions(M, obs_xy, obs_rho)
% Same Hamiltonian as optimizeBVP.m, with one barrier radius per circle. Cached, so an
% angle sweep (same obstacles every solve) builds it once.
persistent key fns
k = [obs_xy(:); obs_rho(:)].';
if ~isempty(fns) && isequal(key, k)
    [ode_rhs_fn, u_star_fn, H_val_fn] = fns{:};
    return
end
fprintf('Building the Hamiltonian symbolically (%d barrier circles)...\n', numel(obs_rho));

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 real
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 real
syms u1 u2 u3 real
syms w_t k_wall k_obs real
x_sym = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];
lam_sym = [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10];
u_sym = [u1; u2; u3];

f_sym = dynamics(x_sym, u_sym);
sigma = jacobian(f_sym, u_sym).' * lam_sym;
a_sym = -sigma .* M.u_max;
u_star = M.u_max .* tanh(a_sym);              % smooth saturating control

x_mid = 0.5*(M.x_min + M.x_max);  x_half = 0.5*(M.x_max - M.x_min);
l_mid = 0.5*(M.l_min + M.l_max);  l_half = 0.5*(M.l_max - M.l_min);
V_lim = k_wall * (((x1 - x_mid)/x_half)^20 + ((x2 - l_mid)/l_half)^20 + ...
                  (x6/M.dx_max)^20 + (x7/M.dl_max)^20 + (x8/M.dth_max)^20);

r_loc = x1 + x2*sin(x4)*cos(x5);  t_loc = x2*sin(x5);
Xp = r_loc*cos(x3) - t_loc*sin(x3);  Yp = r_loc*sin(x3) + t_loc*cos(x3);
V_obs = sym(0);
for j = 1:size(obs_xy, 2)
    V_obs = V_obs + (obs_rho(j)^2 / ((Xp - obs_xy(1, j))^2 + (Yp - obs_xy(2, j))^2))^6;
end
V_obs = k_obs * V_obs;

L_eff = sum(a_sym.*tanh(a_sym) - log(cosh(a_sym)));   % effort cost that u_star minimizes
f_opt = subs(f_sym, u_sym, u_star);
H_sym = w_t + L_eff + V_lim + V_obs + lam_sym.' * f_opt;
ode_rhs = [f_opt; -jacobian(H_sym, x_sym).'];

ode_rhs_fn = matlabFunction(ode_rhs, 'Vars', {x_sym, lam_sym, w_t, k_wall, k_obs});
u_star_fn  = matlabFunction(u_star,  'Vars', {x_sym, lam_sym});
H_val_fn   = matlabFunction(H_sym,   'Vars', {x_sym, lam_sym, w_t, k_wall, k_obs});
key = k;  fns = {ode_rhs_fn, u_star_fn, H_val_fn};
end
