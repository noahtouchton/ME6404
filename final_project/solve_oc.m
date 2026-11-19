%% solve_oc.m - Optimal-control loop around the pin (Pontryagin's principle + bvp4c)
% Optimal control finds the whole trajectory at once: the 11 states, their 11 costates and
% the input U*(t) = u_max*tanh(...) (model/EOM.m) solve a two-point boundary-value problem,
% with the crane at rest at its start position at both ends.
%   1. Fixed time (t_f = 8 s). bvp4c needs a good guess, so start with the pin at 0.80 m,
%      where the trolley alone can carry the payload around it (guess: the payload circling
%      the pin), then move the pin out to 0.90 m a step at a time, each solution being the
%      next guess. Soft walls that are still broken get 3x stiffer.  -> results/oc_fixed.mat
%   2. Minimum time. Free t_f (time scaled to [0, 1], H(t_f) = 0 added) and raise the price
%      of time W_t until t_f stops dropping. J = W_t*t_f + int(L) dt, so as W_t grows the
%      effort term fades and the loop approaches minimum time under the acceleration limits.
%      Saved whenever it beats the fastest loop so far, so Ctrl-C keeps the best.
%                                                                  -> results/oc_min_time.mat
% run_sim.m shows a result; compare_methods.m compares it with the input-shaped loops.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'model', 'generated'));
t_start = tic;

%% Settings
cfg = crane.task();
cfg.r1 = 0.2;  cfg.r2 = 0.5;  cfg.r3 = 0.2;   % effort weights: trolley, hoist, slew
cfg.W_obs = 10;          % pin wall
cfg.k_x   = 1;           % trolley-limit wall
cfg.k_v   = 1;           % velocity walls (must stay >> W_t or speed beats the limits)
cfg.W_t   = 0;           % price of time: 0 while t_f is fixed
t_f       = 8;           % fixed horizon (s). Longer is not easier: the costates get stiff
R_p_start = 0.80;        % first pin position: reach needed 0.85 m < xmax
min_time  = true;        % false: stop after the fixed-time loop
W_t_max   = 1e4;         % stop raising the price of time here
opts  = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'NMax', 20000);
optsT = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'NMax', 50000);

%% 1. Fixed time: continuation in the pin position
R_p_final = cfg.R_p;  cfg.R_p = R_p_start;
T_sw = 2*pi*sqrt(cfg.l0/cfg.g);                                       % swing period (s)
sol  = bvpinit(linspace(0, t_f, round(20*t_f/T_sw)), @(t) loop_guess(t, t_f, cfg));
rng(1);  sol.y(12:22, :) = 1e-3*randn(11, numel(sol.x));   % costates exactly 0 make the Newton matrix singular
[sol, ok, info] = solve_bvp(sol, cfg, false, opts);
fprintf('R_p = %.4f  %s\n', cfg.R_p, info);
if ~ok, error('First solve failed. Try another t_f, a smaller W_obs or a smaller R_p_start.'); end
step = 0.01;
while cfg.R_p < R_p_final - 1e-9
    c = cfg;  c.R_p = min(cfg.R_p + step, R_p_final);
    [s, ok, info] = solve_bvp(sol, c, false, opts);
    fprintf('R_p = %.4f  %s\n', c.R_p, info);
    if ok
        sol = s;  cfg = c;  step = min(1.5*step, 0.02);   % accept: it becomes the next guess
    else
        step = step/2;                                    % retry closer to the last solution
        if step < 1e-4, error('Stalled at R_p = %.4f', cfg.R_p); end
    end
end
[sol, cfg, ~, M] = harden(sol, cfg, false, opts, 6);
S = pack(sol, cfg, false, 'OC fixed 8 s', toc(t_start));
save(crane.file('oc_fixed.mat'), '-struct', 'S');
crane.report('Fixed time', M, cfg);
if ~M.admissible, warning('The fixed-time loop breaks: %s', strjoin(M.why, ', ')); end
if ~min_time, return; end

%% 2. Minimum time: free t_f, raise the price of time
% dJ*/dt_f = H(t_f), so W_t = -H(t_f) is the price at which the current t_f is already optimal
cfg.W_t = max(-hamiltonian_value(deval(sol, t_f), crane.params(cfg)), 1e-4);
[solT, ok, info] = solve_bvp(bvpinit(sol.x/t_f, @(tau) deval(sol, tau*t_f), t_f), cfg, true, optsT);
if ~ok, error('Free-t_f start failed: %s', info); end
M = check(solT, cfg, true);
if ~M.admissible, error('The starting loop breaks: %s', strjoin(M.why, ', ')); end
S = pack(solT, cfg, true, 'OC min-time', toc(t_start));
save(crane.file('oc_min_time.mat'), '-struct', 'S');
fprintf('start: W_t = %8.3g  t_f = %6.3f s\n', cfg.W_t, solT.parameters);
fac = 2;  n_small = 0;  best_tf = solT.parameters;
while cfg.W_t < W_t_max
    c = cfg;  c.W_t = min(cfg.W_t*fac, W_t_max);
    [s, ok, info] = solve_bvp(solT, c, true, optsT);
    if ok, [s, c, ok, M] = harden(s, c, true, optsT, 4); end
    if ok && M.admissible
        gain = 1 - s.parameters/solT.parameters;
        solT = s;  cfg = c;
        fprintf('W_t = %8.3g  t_f = %6.3f s (%+.2f %%)\n', cfg.W_t, solT.parameters, -100*gain);
        if solT.parameters < best_tf                      % stiffer walls can slow it: keep the fastest
            best_tf = solT.parameters;
            S = pack(solT, cfg, true, 'OC min-time', toc(t_start));
            save(crane.file('oc_min_time.mat'), '-struct', 'S');
        end
        if gain < 2e-3, n_small = n_small + 1; else, n_small = 0; end
        if n_small >= 2, fprintf('t_f converged.\n'); break; end   % two steps in a row gained < 0.2 %
        fac = min(fac*1.5, 8);
    else
        if ok, info = strjoin(M.why, ', '); end
        fprintf('   W_t = %.3g rejected (%s)\n', c.W_t, strtrim(info));
        fac = sqrt(fac);                                  % smaller jump next time
        if fac < 1.02, fprintf('Stalled: no larger W_t gives an admissible loop.\n'); break; end
    end
end
S = load(crane.file('oc_min_time.mat'));
crane.report('Minimum time (fastest found)', check(S.sol, S.cfg, true), S.cfg);

%% Local functions
function [sol, ok, info] = solve_bvp(guess, cfg, free, opts)
% One bvp4c solve with every failure check; returns the guess unchanged on failure.
% 22 boundary conditions: at rest at X0 at both ends, eta(0) = 0, and the costate of eta zero
% at the end. eta(t_f) = 2*pi*kdir is not imposed (that condition is degenerate); the loop
% count comes from the guess and the pin wall and is checked here instead.
% Free t_f: time scaled to tau = t/t_f in [0, 1], t_f = sol.parameters, and H(t_f) = 0 added.
    p = crane.params(cfg);  X0 = cfg.X0;
    bc = @(ya, yb) [ya(1:10) - X0; ya(11); yb(1:10) - X0; yb(22)];
    if free
        odefun = @(tau, Z, T) T * hamiltonian_ode(T*tau, Z, p);
        bcfun  = @(ya, yb, T) [bc(ya, yb); hamiltonian_value(yb, p)];
    else
        odefun = @(t, Z) hamiltonian_ode(t, Z, p);
        bcfun  = bc;
    end
    lastwarn('');
    try
        sol = bvp4c(odefun, bcfun, guess, opts);
    catch err                                   % e.g. "singular Jacobian"
        sol = guess;  ok = false;  info = err.message;  return
    end
    msg = lastwarn;                             % non-empty: bvp4c complained
    Zf = sol.y(:, end);
    end_err = max(abs([Zf(1:10) - X0; Zf(11) - 2*pi*cfg.kdir]));
    ok = isempty(msg) && sol.stats.maxres < 1e-3 && end_err < 1e-3 && (~free || sol.parameters > 0);
    info = sprintf('maxres=%.1e end_err=%.1e mesh=%d %s', sol.stats.maxres, end_err, numel(sol.x), msg);
    if free, info = sprintf('t_f=%.3f s %s', sol.parameters, info); end
    if ~ok, sol = guess; end
end

function [sol, cfg, ok, M] = harden(sol, cfg, free, opts, rounds)
% Stiffen (x3) every soft wall the loop breaks and solve again, until every limit holds
    ok = true;
    for n = 0:rounds
        M = check(sol, cfg, free);
        if M.admissible || n == rounds || any(strcmp(M.why, 'windings')), return; end
        c = cfg;
        if any(strcmp(M.why, 'velocity')), c.k_v = 3*c.k_v; end
        if any(strcmp(M.why, 'x')),        c.k_x = 3*c.k_x; end
        if any(strcmp(M.why, 'payload-pin') | strcmp(M.why, 'cable-pin')), c.W_obs = 3*c.W_obs; end
        if isequal(c, cfg), return; end                  % nothing left to stiffen
        fprintf('   stiffen (%s): k_v %.3g  k_x %.3g  W_obs %.3g\n', strjoin(M.why, ', '), c.k_v, c.k_x, c.W_obs);
        [s, ok] = solve_bvp(sol, c, free, opts);
        if ~ok, return; end
        sol = s;  cfg = c;
    end
end

function M = check(sol, cfg, free)
% Every limit and clearance of a bvp4c solution, densely sampled in real time
    S = pack(sol, cfg, free, '', 0);
    M = crane.metrics(crane.trajectory(S, linspace(0, S.t_f, 1500)), cfg);
end

function S = pack(sol, cfg, free, label, compute_time)
% A saved solution (see crane.m)
    if free, t_f = sol.parameters; else, t_f = sol.x(end); end
    S = struct('cfg', cfg, 'kind', 'optimal', 'sol', sol, 'time_scaled', free, 't_f', t_f, 'label', label, ...
               'source', 'solve_oc', 'created', char(datetime('now')), 'compute_time', compute_time);
end

function Y = loop_guess(t, t_f, cfg)
% Initial guess: the payload circles the pin |kdir| times hanging straight down, starting and
% stopping at rest. Costates are filled in by the caller.
    tau = t / t_f;
    s   = tau - sin(2*pi*tau)/(2*pi);             % 0 -> 1 with zero slope at both ends
    sd  = (1 - cos(2*pi*tau)) / t_f;              % ds/dt
    pin = [cfg.R_p; 0];
    rho = norm([cfg.x0; 0] - pin);                % loop radius: the start's distance from the pin
    a   = 2*pi*abs(cfg.kdir)*s;                   % angle around the pin
    Q   = pin + rho*[-cos(a); -sign(cfg.kdir)*sin(a)];
    Qd  = rho*2*pi*abs(cfg.kdir)*sd*[sin(a); -sign(cfg.kdir)*cos(a)];
    x   = norm(Q);  th = atan2(Q(2), Q(1));       % trolley right above the payload
    dx  = dot(Q, Qd) / x;
    dth = (Q(1)*Qd(2) - Q(2)*Qd(1)) / x^2;
    Y   = [x; cfg.l0; th; 0; 0; dx; 0; dth; 0; 0; cfg.kdir*2*pi*s; zeros(11, 1)];
end
