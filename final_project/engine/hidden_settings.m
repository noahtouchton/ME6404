function H = hidden_settings()
% Everything the lab doesn't care about: hardware limits, cost weights, wall stiffness,
% solver tolerances and continuation schedules. lab_tune.m never needs to touch these.
% The solver hardens walls by itself (and saves the hardened values), so these are only
% starting values.

% --- Hardware limits. Fixed. ---
H.g       = 9.81;
H.xmin    = 0;       H.xmax   = 0.925;     % trolley radius (m)
H.lmin    = 0;       H.lmax   = 1.450;     % cable length (m)
H.dxmax   = 72.6e-3; H.dlmax  = 30e-3;     % m/s
H.dthmax  = deg2rad(32.4);                 % rad/s
H.ddxmax  = 72.6e-3; H.ddlmax = 30e-3;     % m/s^2
H.ddthmax = deg2rad(32.4);                 % rad/s^2
H.drive_tau = 0.05;                        % s, velocity-loop time constant used when replaying the upload file

% --- Cost weights (relative effort of the three drives) ---
H.r1 = 0.2;   H.r2 = 0.5;   H.r3 = 0.2;    % trolley, hoist, slew: hoist is the "expensive" one
H.W_obs = 10;                              % pin wall
H.k_x   = 1;                               % position walls
H.k_v   = 1;                               % velocity walls

% --- Admissibility tolerances (soft walls may graze a limit by this fraction) ---
H.tol = 1e-3;

% --- bvp4c ---
% Loose on purpose: what's uploaded is checked by simulate_commands, which integrates the
% real dynamics with the exact file commands. At 1e-2 that replay is as clean as at 1e-3
% (same pin gap, 0.3 mm residual, 0.05 mm end error), and a re-solve takes ~15 s instead of ~4 min.
H.bvp = bvpset('RelTol', 1e-2, 'AbsTol', 1e-3, 'NMax', 15000);
H.maxres = 1e-2;              % accept a solve only if its residual is below this
H.warm_mesh = 600;            % resample a warm start onto this many points (bvp4c adds where needed)

% --- Continuation (homotopy from the last good solution to the requested problem) ---
H.step0     = 0.25;           % first fraction of the way to the target
H.step_min  = 1/512;          % give up below this
H.max_hard  = 4;              % wall-hardening rounds per step
H.hard_fac  = 3;              % wall weight multiplier per round

% --- Minimum time (W_t continuation: min t_f + (1/W_t)∫L dt -> pure min time as W_t -> ∞) ---
H.W_t_max   = 1e4;
H.W_t_fac0  = 2;
H.tf_tol    = 2e-3;           % converged after two accepted steps that each cut t_f < 0.2 %

% --- Cold start ---
H.guess_margin = 0.03;        % guess path passes this far outside the clearance circle (m)
H.guess_slack  = 1.4;         % guess runs each axis at most 1/1.4 of its velocity limit
H.seed = 1;                   % rng seed for the tiny random costate guess (repeatable runs)
end
