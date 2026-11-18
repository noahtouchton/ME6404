function R = simulate_commands(S, cmd, t_watch)
% Open-loop check of the exported file: drive the crane model with the velocity commands
% exactly as uploaded (each drive tracks its command with time constant drive_tau, never
% exceeding the acceleration limits), then keep watching t_watch seconds after the stop.
% This integrates the real dynamics independently of the BVP solver, so it's the honest
% test of what goes to the crane.
if nargin < 3, t_watch = 5; end
cfg = S.cfg;  hw = cfg.hw;  p = crane_params(cfg);
dt  = 0.005;
t_cmd = cmd.t(end);
t   = (0:dt:t_cmd + t_watch).';
N   = numel(t);
amax = [hw.ddxmax; hw.ddlmax; hw.ddthmax];             % model input order: [ẍ; l̈; θ̈]

X = zeros(11, N);  X(1:10, 1) = cfg.X_start;
U = zeros(3, N);
for k = 1:N-1
    j = find(cmd.t <= t(k) + 1e-9, 1, 'last');          % active row of the file
    v_cmd = cmd.v(j, [1 3 2]).' * (t(k) < t_cmd);       % file order [x θ l] -> model order [x l θ]
    u_of  = @(x) max(min((v_cmd - x(6:8)) / cfg.drive_tau, amax), -amax);
    U(:, k) = u_of(X(:, k));
    f  = @(x) crane_dynamics(x, u_of(x), p);           % RK4 with the drive law inside
    k1 = f(X(:, k));            k2 = f(X(:, k) + dt/2*k1);
    k3 = f(X(:, k) + dt/2*k2);  k4 = f(X(:, k) + dt*k3);
    X(:, k+1) = X(:, k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

R.t = t;  R.X = X;  R.U = U;  R.t_cmd = t_cmd;
R.G = crane_geometry(X, cfg);
Xt  = cfg.X_end;  Xt(11) = 0;
Qt  = crane_geometry(Xt, cfg).Qp(1:2);                  % where the payload should hang at the end
after = t >= t_cmd;
R.offset_mm       = 1000 * vecnorm(R.G.Qp(1:2, :) - Qt);   % payload distance from its rest point
R.residual_peak_mm = max(R.offset_mm(after));            % residual swing after the stop
R.end_pose_err    = [1000*(X(1, end) - Xt(1)), rad2deg(X(3, end) - Xt(3)), 1000*(X(2, end) - Xt(2))];  % mm, deg, mm
R.d_pin_min_mm    = 1000 * min(R.G.d_pin);
R.d_cable_min_mm  = 1000 * min(R.G.d_cable);
end
