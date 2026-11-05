%% find_min_time.m - Shrink t_f to the fastest admissible loop around the pin
% Reads crane_solution.mat, frees the final time and raises the time price W_t
% until t_f stops dropping. After every accepted step the file is overwritten with
% the fastest admissible solution so far (Ctrl-C keeps the best result).
%
% Theory: J = W_t*t_f + int(L) dt  <=>  min t_f + (1/W_t)*int(L) dt, so as
% W_t -> inf the effort term vanishes and the solution approaches minimum time
% under the acceleration limits (bang-bang in the limit). The soft walls (pin,
% velocities, trolley) must keep costing much more than W_t, so any wall that
% gets violated is hardened (x3) and the step is re-solved.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'), fullfile(here, 'model', 'generated'));
io_file  = fullfile(here, 'results', 'crane_solution.mat');
bak_file = fullfile(here, 'results', 'crane_solution_backup.mat');

W_t_max  = 1e4;     % stop here even if t_f is still dropping
tf_tol   = 2e-3;    % converged once two accepted steps in a row each cut t_f by < 0.2 %
max_hard = 4;       % wall-hardening attempts per W_t step
optsT    = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'NMax', 50000);

%% 1. Load and set up the free-tf problem
if ~isfile(io_file), error('No %s yet: run find_path.m first.', io_file); end
S   = load(io_file);
cfg = S.cfg;
copyfile(io_file, bak_file);                 % undo: rename the backup back
fprintf('Loaded %s (%s, t_f = %.2f s). Backup: %s\n', io_file, S.source, S.t_f, bak_file);

if S.time_scaled
    solT = S.sol;                            % already free-tf: resume from it
else
    % Map the fixed-tf solution onto tau in [0,1]. dJ*/dt_f = H(t_f), so W_t = -H
    % is exactly the time price at which the current t_f is already optimal.
    H_fix   = hamiltonian_value(deval(S.sol, S.t_f), crane_params(cfg));
    cfg.W_t = max(-H_fix, 1e-4);
    guess   = bvpinit(S.sol.x / S.t_f, @(tau) deval(S.sol, tau*S.t_f), S.t_f);
    [solT, ok, info] = try_solve_T(guess, crane_params(cfg), cfg.X0, cfg.kdir, optsT);
    if ~ok, error('Free-tf start failed: %s', info); end
end
M = metrics_of(solT, cfg);
fprintf('start: W_t = %8.3g  t_f = %6.3f s  %s\n', cfg.W_t, solT.parameters, status_of(M));
if ~M.admissible, error('Starting solution is not admissible (%s); fix find_path first.', strjoin(M.why, ', ')); end

%% 2. Raise W_t until t_f converges
history = [cfg.W_t, solT.parameters];        % [W_t, t_f] of every accepted step
fac = 2;  n_small = 0;
while cfg.W_t < W_t_max
    c = cfg;  c.W_t = min(cfg.W_t*fac, W_t_max);
    [s_try, ok, info] = try_solve_T(solT, crane_params(c), c.X0, c.kdir, optsT);

    % Re-solve with harder walls until the result respects every limit
    M = struct('admissible', false, 'why', {{'solve failed'}});
    n_hard = 0;
    while ok
        M = metrics_of(s_try, c);
        if M.admissible || n_hard >= max_hard || any(strcmp(M.why, 'windings')), break; end
        if any(strcmp(M.why, 'velocity')),                          c.k_v   = 3*c.k_v;   end
        if any(strcmp(M.why, 'x')),                                 c.k_x   = 3*c.k_x;   end
        if any(strcmp(M.why, 'payload-pin') | strcmp(M.why, 'cable-pin')), c.W_obs = 3*c.W_obs; end
        n_hard = n_hard + 1;
        fprintf('   harden (%s): k_v = %.3g  k_x = %.3g  W_obs = %.3g\n', strjoin(M.why, ', '), c.k_v, c.k_x, c.W_obs);
        [s_try, ok, info] = try_solve_T(s_try, crane_params(c), c.X0, c.kdir, optsT);
    end

    if ok && M.admissible
        gain = 1 - s_try.parameters / solT.parameters;
        solT = s_try;  cfg = c;
        history(end+1, :) = [cfg.W_t, solT.parameters];      %#ok<SAGROW>
        fprintf('W_t = %8.3g  t_f = %6.3f s (-%5.2f %%)  %s\n', cfg.W_t, solT.parameters, 100*gain, status_of(M));
        save_result(io_file, solT, cfg, history);
        if gain < tf_tol, n_small = n_small + 1; else, n_small = 0; end
        if n_small >= 2, fprintf('t_f converged.\n'); break; end
        fac = min(fac*1.5, 8);
    else
        if ~ok, reason = info; else, reason = strjoin(M.why, ', '); end
        fprintf('   W_t = %.3g rejected (%s)\n', c.W_t, strtrim(reason));
        fac = sqrt(fac);                                      % smaller jump next time
        if fac < 1.02, fprintf('Stalled: no larger W_t gives an admissible solution.\n'); break; end
    end
end

%% 3. Summary
fprintf('\nFastest admissible maneuver: t_f = %.3f s  (W_t = %.3g, k_v = %.3g, k_x = %.3g, W_obs = %.3g)\n', ...
    solT.parameters, cfg.W_t, cfg.k_v, cfg.k_x, cfg.W_obs);
fprintf('Saved to %s. Run run_sim.m to see it.\n', io_file);

%% Local functions
function M = metrics_of(solT, cfg)
% Admissibility checks for a free-tf solution (densely sampled in real time).
    S = struct('cfg', cfg, 'sol', solT, 'time_scaled', true, 't_f', solT.parameters);
    M = crane_metrics(crane_trajectory(S, linspace(0, solT.parameters, 1500)), cfg);
end

function s = status_of(M)
% One-line summary of how hard the limits are being used.
    s = sprintf('|v|/vmax = %.2f %.2f %.2f  |U|/Umax = %.2f %.2f %.2f  pin %.1f mm  cable %.1f mm', ...
        M.v_ratio, M.u_ratio, 1000*M.d_pin_min, 1000*M.d_cable_min);
end

function save_result(io_file, solT, cfg, history)
% Overwrite the solution file with a free-tf result (run_sim reads the same fields).
    S = struct('cfg', cfg, 'sol', solT, 'time_scaled', true, 't_f', solT.parameters, ...
               'source', 'find_min_time', 'created', char(datetime('now')), 'W_t_history', history);
    save(io_file, '-struct', 'S');
end
