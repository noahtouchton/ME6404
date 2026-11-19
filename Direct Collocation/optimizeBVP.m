clc; clear; close all;

%% Limits
x_min = 0.0; x_max = 0.925;
l_min = 0.0; l_max = 1.450;
dx_min = -72.6e-3; dx_max = -dx_min;
dl_min = -30.0e-3; dl_max = -dl_min;
dth_min = -deg2rad(32.4); dth_max = -dth_min;
ddx_min = -72.6e-3; ddx_max = -ddx_min;
ddth_min = -deg2rad(32.4); ddth_max = -ddth_min;
ddl_min = -30.0e-3; ddl_max = -ddl_min;

u_max = [ddx_max; ddl_max; ddth_max];
R_diag = 1 ./ (u_max.^2); % Bryson's Rule (the effort cost below reduces to 0.5*u'*R*u for small u)

%% Starting and Ending States

% State = [x; l; th; phi; psi; dx; dl; dth; dphi; dpsi]

% Trolley
x_start = 0.8; x_end = 0.4;
dx_start = 0.0; dx_end = 0.0;

% Hoist
l_start = 1.0; l_end = 1.2;
dl_start = 0.0; dl_end = 0.0;

% Slew
th_start = deg2rad(-20.0); th_end = deg2rad(30.0);
dth_start = 0.0; dth_end = 0.0;

% Phi and Psi
phi_start = 0.0; phi_end = 0.0;
psi_start = 0.0; psi_end = 0.0;
dphi_start = 0.0; dphi_end = 0.0;
dpsi_start = 0.0; dpsi_end = 0.0;

xInitial = [x_start; l_start; th_start; phi_start; psi_start; dx_start; dl_start; dth_start; dphi_start; dpsi_start];
xFinal = [x_end; l_end; th_end; phi_end; psi_end; dx_end; dl_end; dth_end; dphi_end; dpsi_end];

%% Obstacle and Geometry Parameters

H_jib = 2.0;
R_pay = 0.08;

% The real obstacle: the same wall as optimize.m (used by sim_trajectory.m to draw and check it)
obs(1).type   = 'box';
obs(1).center = [0.3; 0.05; 0.750];
obs(1).half   = [0.35; 0.02; 0.750];

% PMP needs a smooth barrier, so the wall is modeled as a chain of vertical circles
% along its centerline (radius = half thickness + payload radius). Z is ignored: the
% payload hangs at Z = H_jib - l <= 1.0 m here, below the 1.5 m wall top.
obs_rho = obs(1).half(2) + R_pay;
x_ends  = obs(1).center(1) + [-1 1]*(obs(1).half(1) - obs(1).half(2));
n_circ  = ceil(diff(x_ends) / 0.05) + 1;                          % <= 5 cm apart (~3 mm sag between circles)
obs_xy  = [linspace(x_ends(1), x_ends(2), n_circ); obs(1).center(2)*ones(1, n_circ)];

check_endpoint_clearance(xInitial, xFinal, obs, R_pay, H_jib);

%% Symbolically Generate Costate and Hamiltonian Functions
[ode_rhs_fn, u_star_fn, H_val_fn] = build_pmp_symbols( ...
    x_min, x_max, l_min, l_max, dx_max, dl_max, dth_max, ...
    u_max, obs_xy, obs_rho);

%% Build Initial Guess Mesh (waypoints around the wall's outer end, zero costates)
N_init = 60;
tau_init = linspace(0, 1, N_init);
tf_guess = 12.0;

% Payload XY points the guess passes through (payload hanging straight down)
wp_XY = [0.82, 0.75;        % 1: past the wall's outer end (0.65 + rho + margin)
         0.05, 0.22];       % 2: up and away from the wall before heading in
x_guess = waypoint_guess(xInitial, xFinal, wp_XY, tau_init, tf_guess);

lam_guess = zeros(10, N_init);
Y_guess = [x_guess; lam_guess];

bvp_opts = bvpset('RelTol', 1e-3, 'AbsTol', 1e-5, 'NMax', 15000, 'Stats', 'off');

%% Stage 0: Fixed-tf Cold Start (Finds true costates lambda(t) first)
% The obstacle is on from the start: the guess is already around the wall, and the
% barrier keeps the solution there instead of letting it cut through.
fprintf('Solving PMP Stage 0 (Fixed tf = %.1f s, finding costates)...\n', tf_guess);
p_stage0 = [0.20, 1e-6, 1e-2];
sol0 = bvpinit(tau_init, @(tau) interp1(tau_init, Y_guess', tau)');
lastwarn('');
sol0 = bvp4c(@(tau, Y) bvp_ode(tau, Y, tf_guess, ode_rhs_fn, p_stage0), ...
             @(Ya, Yb) [Ya(1:10) - xInitial; Yb(1:10) - xFinal], ...
             sol0, bvp_opts);
if ~isempty(lastwarn), error('Stage 0 did not converge: %s', lastwarn); end

sol = bvpinit(sol0.x, @(tau) deval(sol0, tau), tf_guess);

%% Continuation Loop

% Each stage is reached in adaptive sub-steps from the previous one. A sub-step is
% accepted only if bvp4c met its tolerance; otherwise the step is halved. This stops a
% failed solve from being handed to the next stage as its warm start.

% Columns: [w_t, k_wall, k_obs]
stages = [
    0.20, 1e-5, 1e-2;   % Stage 1: Free tf, gentle speed
    0.50, 1e-4, 1e-1;   % Stage 2: Firmer walls and obstacle
    1.50, 5e-4, 1.0;    % Stage 3: Obstacle effectively solid & move faster
    3.00, 1e-3, 1.0;    % Stage 4: Final fast, near-minimum-time run
    3.00, 1.0,  1.0     % Stage 5: Harden the limit walls (at k_wall = 1e-3 the power-20 wall
                        %          still lets a velocity reach ~140% of its limit)
];

p_prev = p_stage0;
for s = 1:size(stages, 1)
    p_to = stages(s, :);
    fprintf('Solving PMP Stage %d/%d (w_t = %.2f, k_wall = %.1e, k_obs = %.1e)...\n', ...
            s, size(stages, 1), p_to(1), p_to(2), p_to(3));
    a = 0;  step = 1;                         % try the whole stage at once first
    while a < 1
        a_try = min(1, a + step);
        p_try = blend(p_prev, p_to, a_try);
        lastwarn('');
        try
            sol_try = bvp4c(@(tau, Y, tf) bvp_ode(tau, Y, tf, ode_rhs_fn, p_try), ...
                @(Ya, Yb, tf) bvp_bc(Ya, Yb, tf, xInitial, xFinal, H_val_fn, p_try), ...
                sol, bvp_opts);
            ok = isempty(lastwarn);
        catch err
            ok = false;  fprintf('    (%s)\n', err.message);
        end
        if ok
            sol = sol_try;  a = a_try;  step = min(2*step, 1);
            if a < 1, fprintf('    %3.0f%% of the way, tf = %.2f s\n', 100*a, sol.parameters); end
        else
            step = step/2;
            fprintf('    step failed, halving to %.3f\n', step);
            if step < 1/64, error('Stuck in stage %d at %.0f%% of the way', s, 100*a); end
        end
    end
    p_prev = p_to;
    fprintf('  Stage %d done: tf = %.2f s, %d mesh points\n', s, sol.parameters, numel(sol.x));
end

%% Unpack Adaptive Mesh Solution and Save trajectory

tf_opt = sol.parameters;
t_grid = sol.x * tf_opt;
x_opt = sol.y(1:10, :);
lam_opt = sol.y(11:20, :);

u_opt = zeros(3, length(t_grid));
for k=1:length(t_grid)
    u_opt(:,k) = u_star_fn(x_opt(:,k), lam_opt(:,k));
end

H_t = arrayfun(@(k) H_val_fn(x_opt(:,k), lam_opt(:,k), p_prev(1), p_prev(2), p_prev(3)), 1:length(t_grid));
fprintf('Converged! Optimal tf = %.4f s (Adaptive mesh nodes = %d)\n', tf_opt, length(t_grid));
fprintf('  PMP check: H(t) spread %.1e (should be ~0 for free tf)\n', max(abs(H_t)));
fprintf('  peak velocity / limit: trolley %.0f%%, hoist %.0f%%, slew %.0f%%\n', ...
        100*max(abs(x_opt(6,:)))/dx_max, 100*max(abs(x_opt(7,:)))/dl_max, 100*max(abs(x_opt(8,:)))/dth_max);

x_lb = [x_min; l_min; -inf; -inf; -inf; dx_min; dl_min; dth_min; -inf; -inf];
x_ub = [x_max; l_max; inf; inf; inf; dx_max; dl_max; dth_max; inf; inf];

traj.t       = t_grid;
traj.x       = x_opt;
traj.u       = u_opt;
traj.lam     = lam_opt;
traj.tf      = tf_opt;
traj.w_t     = stages(end, 1);
traj.N       = length(t_grid) - 1;
traj.x_lb    = x_lb;    traj.x_ub = x_ub;      % limits, so the sim can check them
traj.u_lb    = -u_max;  traj.u_ub = u_max;
traj.obs     = obs;                            % the real wall, drawn and checked by the sim
traj.H_jib   = H_jib;
traj.R_pay   = R_pay;
traj.created = char(datetime('now'));
out_file = fullfile(fileparts(mfilename('fullpath')), 'trajectory.mat');
save(out_file, '-struct', 'traj');
fprintf('Saved trajectory to %s\n', out_file);

% Inputs as % of max acceleration, one row per mesh point (overwritten every run)
export_xlsx(fullfile(fileparts(mfilename('fullpath')), 'inputs_bvp.xlsx'), t_grid, u_opt, u_max);

%% BVP Functions

function dYd_tau = bvp_ode(tau, Y, tf, ode_rhs_fn, p_stage)
    x = Y(1:10);
    lam = Y(11:20);
    dYd_tau = tf * ode_rhs_fn(x, lam, p_stage(1), p_stage(2), p_stage(3));
end

function res = bvp_bc(Ya, Yb, ~, xInitial, xFinal, H_val_fn, p_stage)
% 21 Boundary Conditions: 10 initial states + 10 final states + H(tf) = 0
xa = Ya(1:10);
xb = Yb(1:10);
lamb = Yb(11:20);

H_tf = H_val_fn(xb, lamb, p_stage(1), p_stage(2), p_stage(3));
res = [xa - xInitial; ...
    xb - xFinal; ...
    H_tf];
end

function p = blend(p_from, p_to, a)
% Continuation between two stages: w_t linearly, the stiffness weights geometrically
% (they span orders of magnitude, so a linear blend would jump almost all the way at once).
p = p_from + a*(p_to - p_from);
for i = 2:3
    if p_from(i) > 0 && p_to(i) > 0
        p(i) = p_from(i) * (p_to(i)/p_from(i))^a;
    end
end
end

function x_guess = waypoint_guess(xInitial, xFinal, wp_XY, tau, tf)
% State guess: payload XY goes start -> waypoints -> finish along straight segments,
% with a smoothstep time profile (rest at both ends); cable length blends linearly in s.
% Trolley radius/slew come from the payload XY (payload hanging straight down).
s  = 3*tau.^2 - 2*tau.^3;
P0 = xInitial(1) * [cos(xInitial(3)); sin(xInitial(3))];
P1 = xFinal(1)   * [cos(xFinal(3));   sin(xFinal(3))];
P  = [P0, wp_XY, P1];
cum = [0, cumsum(vecnorm(diff(P, 1, 2)))];  cum = cum / cum(end);
XY = interp1(cum, P.', s).';

x_guess = zeros(10, numel(tau));
x_guess(1, :) = vecnorm(XY);
x_guess(3, :) = unwrap(atan2(XY(2, :), XY(1, :)));
x_guess(3, :) = x_guess(3, :) - x_guess(3, 1) + xInitial(3);
x_guess(2, :) = xInitial(2) + (xFinal(2) - xInitial(2)) * s;
for i = 1:3
    x_guess(i+5, :) = gradient(x_guess(i, :), tau) / tf;         % d/dt = (1/tf) d/dtau
end
end

function check_endpoint_clearance(xInitial, xFinal, obs, R_pay, H_jib)
% The barrier is infinite inside the obstacle, so both rest poses must be outside it.
names = {'start', 'end'};  X = [xInitial, xFinal];
for k = 1:2
    p = [X(1,k)*cos(X(3,k)); X(1,k)*sin(X(3,k)); H_jib - X(2,k)];
    for j = 1:numel(obs)
        d = obstacle_sdf(p, obs(j), R_pay);
        if d <= 0
            error('The %s pose is inside obstacle %d (clearance %.1f mm). Move the pose or the obstacle.', ...
                  names{k}, j, 1000*d);
        end
    end
end
end

%% ========================================================================
%  SYMBOLIC HAMILTONIAN & COSTATE GENERATOR
%  ========================================================================
function [ode_rhs_fn, u_star_fn, H_val_fn] = build_pmp_symbols( ...
    x_min, x_max, l_min, l_max, dx_max, dl_max, dth_max, ...
    u_max, obs_xy, obs_rho)

    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 real
    syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 real
    syms u1 u2 u3 real
    syms w_t k_wall k_obs real

    x_sym   = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];
    lam_sym = [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10];
    u_sym   = [u1; u2; u3];

    % 1. Call your existing dynamics(x, u) symbolically
    f_sym = dynamics(x_sym, u_sym);

    % 2. Extract control matrix B(x) = df/du and switching function sigma = B'*lam
    B_sym     = jacobian(f_sym, u_sym);
    sigma_sym = B_sym.' * lam_sym;

    % 3. Smooth saturating control u*(x, lam) = u_max * tanh(-sigma * u_max)
    % It minimizes H only for the effort cost in step 6, which is why that cost is not
    % 0.5*u'*R*u: with the quadratic cost dH/du ~= 0 at this u* and the costate
    % equations below would not be the PMP conditions of any problem.
    a_sym = -sigma_sym .* u_max;
    u_star_sym = u_max .* tanh(a_sym);

    % 4. Smooth polynomial wall penalties (power 20) for state limits
    x_mid  = 0.5 * (x_min + x_max);  x_half = 0.5 * (x_max - x_min);
    l_mid  = 0.5 * (l_min + l_max);  l_half = 0.5 * (l_max - l_min);

    V_lim = k_wall * ( ((x1 - x_mid)/x_half)^20 + ...
        ((x2 - l_mid)/l_half)^20 + ...
        (x6 / dx_max)^20 + ...
        (x7 / dl_max)^20 + ...
        (x8 / dth_max)^20 );

    % 5. One-sided obstacle barrier: (rho/d)^12 per circle of the wall chain.
    % = k_obs on a circle's surface, -> 0 far away, -> infinity inside.
    r_loc = x1 + x2 * sin(x4) * cos(x5);
    t_loc = x2 * sin(x5);
    Xp    = r_loc * cos(x3) - t_loc * sin(x3);
    Yp    = r_loc * sin(x3) + t_loc * cos(x3);

    V_obs = sym(0);
    for j = 1:size(obs_xy, 2)
        d_sq  = (Xp - obs_xy(1, j))^2 + (Yp - obs_xy(2, j))^2;
        V_obs = V_obs + (obs_rho^2 / d_sq)^6;
    end
    V_obs = k_obs * V_obs;

    % 6. Effort cost that u* actually minimizes: L_i(u) = (1/u_max)[u*atanh(u/u_max) + (u_max/2)*ln(1-(u/u_max)^2)].
    % At u = u_max*tanh(a) this is a*tanh(a) - ln(cosh(a)), and for small u it is 0.5*u^2/u_max^2 (Bryson).
    L_eff     = sum(a_sym .* tanh(a_sym) - log(cosh(a_sym)));
    f_opt_sym = subs(f_sym, u_sym, u_star_sym);
    L_sym     = w_t + L_eff + V_lim + V_obs;
    H_sym     = L_sym + lam_sym.' * f_opt_sym;

    % 7. Derive costate ODEs: dlam/dt = -(dH/dx)^T
    dlam_sym    = -jacobian(H_sym, x_sym).';
    ode_rhs_sym = [f_opt_sym; dlam_sym];

    % Convert symbolic expressions into fast anonymous MATLAB functions
    ode_rhs_fn = matlabFunction(ode_rhs_sym, 'Vars', {x_sym, lam_sym, w_t, k_wall, k_obs});
    u_star_fn  = matlabFunction(u_star_sym,  'Vars', {x_sym, lam_sym});
    H_val_fn   = matlabFunction(H_sym,       'Vars', {x_sym, lam_sym, w_t, k_wall, k_obs});
end
