%% find_path.m - Solve the optimal loop around the pin and save it
% Solves the fixed-horizon TPBVP by continuation and writes crane_solution.mat.
%   find_min_time.m  -> shrinks t_f to the fastest admissible maneuver (updates the file)
%   run_sim.m        -> plots and animates whatever is in the file (no solving)
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'), fullfile(here, 'model', 'generated'));
out_file = fullfile(here, 'results', 'crane_solution.mat');

%% 1. Settings (all saved with the solution, so run_sim / find_min_time use the same ones)
cfg.g = 9.81;

% Hardware limits (fixed)
cfg.xmax    = 0.925;            cfg.lmax   = 1.450;          % m
cfg.dxmax   = 72.6e-3;          cfg.dlmax  = 30.0e-3;        % m/s
cfg.dthmax  = deg2rad(32.4);                                 % rad/s
cfg.ddxmax  = 72.6e-3;          cfg.ddlmax = 30.0e-3;        % m/s^2
cfg.ddthmax = deg2rad(32.4);                                 % rad/s^2

% Pin (fixed in the ground frame on the theta = 0 line)
cfg.r_p   = 0.03;     % pin radius (m)
cfg.c_clr = 0.02;     % clearance: payload size + safety margin (m)
cfg.z_top = -0.75;    % pin top relative to the jib (m). At -0.60 the cable clipped the pin top by 0.1 mm

% Cost weights
cfg.r1 = 0.2;  cfg.r2 = 0.5;  cfg.r3 = 0.2;   % effort: trolley, hoist, slew
cfg.W_obs = 10;    % pin wall
cfg.k_x   = 1;     % trolley-limit wall (hardened below if needed)
cfg.k_v   = 1;     % velocity walls (find_min_time hardens these as it speeds up)
cfg.W_t   = 0;     % time price: 0 here because t_f is fixed

% Maneuver: start = end, hanging at rest inboard of the pin
cfg.x0 = 0.70;  cfg.l0 = 1.00;
cfg.X0 = [cfg.x0; cfg.l0; zeros(8, 1)];
cfg.kdir = 1;          % signed loops around the pin: +1 once CCW, -1 once CW, +2 twice CCW
t_f = 8;               % fixed horizon (s). Longer is NOT easier: the costates get stiff

R_p_start = 0.80;      % reach needed 0.85 < xmax: trolley can do it alone (easy first solve)
R_p_final = 0.90;      % reach needed 0.95 > xmax: payload must swing

%% 2. Initial guess: loop(s) around the pin with the payload hanging straight down
opts   = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'NMax', 20000);
T_sw   = 2*pi*sqrt(cfg.l0/cfg.g);                       % swing period (s)
t_mesh = linspace(0, t_f, round(20*t_f/T_sw));          % ~20 mesh points per swing
cfg.R_p = R_p_start;
sol = bvpinit(t_mesh, @(t) loop_guess(t, t_f, cfg.X0, [cfg.R_p; 0], cfg.kdir));
% Costates exactly 0 make every lambda*df/dX term vanish -> singular Newton matrix.
% A tiny random costate guess breaks that degeneracy (rng fixed so runs repeat).
rng(1);
sol.y(12:22, :) = 1e-3 * randn(11, numel(sol.x));

%% 3. Continuation
% First solve at the easiest pin location
[sol, ok, info] = try_solve(sol, crane_params(cfg), cfg.X0, cfg.kdir, opts);
fprintf('R_p = %.4f  %s\n', cfg.R_p, info);
if ~ok, error('First solve failed. Try a different t_f, smaller W_obs, or smaller R_p_start.'); end

% Move the pin outward until swinging is required
step = 0.01;
while cfg.R_p < R_p_final - 1e-9
    c = cfg;  c.R_p = min(cfg.R_p + step, R_p_final);
    [sol_try, ok, info] = try_solve(sol, crane_params(c), c.X0, c.kdir, opts);
    fprintf('R_p = %.4f  %s\n', c.R_p, info);
    if ok
        sol = sol_try;  cfg = c;            % accept: it becomes the next guess
        step = min(1.5*step, 0.02);
    else
        step = step/2;                      % failed: retry from the last good sol with a smaller step
        if step < 1e-4, error('Stalled at R_p = %.4f', cfg.R_p); end
    end
end

% Harden the trolley wall until the trolley stays inside xmax
while max(deval(sol, linspace(0, t_f, 2000), 1)) > cfg.xmax
    c = cfg;  c.k_x = 3*cfg.k_x;
    [sol_try, ok, info] = try_solve(sol, crane_params(c), c.X0, c.kdir, opts);
    fprintf('k_x = %.0f  %s\n', c.k_x, info);
    if ~ok, error('Could not harden the x wall.'); end
    sol = sol_try;  cfg = c;
end

%% 4. Check and save
S = struct('cfg', cfg, 'sol', sol, 'time_scaled', false, 't_f', t_f, ...
           'source', 'find_path', 'created', char(datetime('now')));
M = crane_metrics(crane_trajectory(S, linspace(0, t_f, 2000)), cfg);
fprintf('\nt_f = %.2f s | max|v|/vmax = %.3f %.3f %.3f | payload-pin %.1f mm | cable-pin %.1f mm | reach over xmax %+.1f mm\n', ...
    t_f, M.v_ratio, 1000*M.d_pin_min, 1000*M.d_cable_min, 1000*M.reach_over);
if ~M.admissible
    warning('Solution violates: %s (saved anyway; see run_sim for details)', strjoin(M.why, ', '));
end
save(out_file, '-struct', 'S');
fprintf('Saved %s\n', out_file);
