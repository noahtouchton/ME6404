    clc; clear;

% x = [x; l; th; phi; psi; dx; dl; dth; dphi; dpsi]
% u = [ddx; ddl; ddth]



N = 30;
tf_guess = 5.0; 
tf_min = 0.5;
tf_max = 30.0;
w_t = 0.05;

%% Limits
x_min = 0.0; x_max = 0.925;
l_min = 0.0; l_max = 1.450;
dx_min = -72.6e-3; dx_max = -dx_min;
dl_min = -30.0e-3; dl_max = -dl_min;
dth_min = -deg2rad(32.4); dth_max = -dth_min;
ddx_min = -72.6e-3; ddx_max = -ddx_min;
ddth_min = -deg2rad(32.4); ddth_max = -ddth_min;
ddl_min = -30.0e-3; ddl_max = -ddl_min;

% no limits on theta, or psi, phi (inlcuding derivatives of these 2)

%% Starting and Ending states

% Trolley
x_start = 0.8; x_end = 0.5;
dx_start = 0.0; dx_end = 0.0;

% Hoist
l_start = 1.0; l_end = 1.2;
dl_start = 0.0; dl_end = 0.0;

% Slew
th_start = deg2rad(-10.0); th_end = deg2rad(30.0);
dth_start = 0.0; dth_end = 0.0;

% Phi and Psi
phi_start = 0.0; phi_end = 0.0;
psi_start = 0.0; psi_end = 0.0;
dphi_start = 0.0; dphi_end = 0.0;
dpsi_start = 0.0; dpsi_end = 0.0;



%% Obstacles (ground frame, Z up from the floor)
H_jib = 2.0;   % height of the jib above the floor
R_pay = 0.08;  % payload radius (clearance margin)

% Obstacle 1: wall/block, 0.4m wide, 0.04m thick, 1.0m tall
obs(1).type   = 'box';
obs(1).center = [0.3; 0.05; 0.750];
obs(1).half   = [0.35; 0.02; 0.750];

% Obstacle 2: pillar at X=0.60, Y=0.10, radius 0.05m, height 0 to 1.2m
%obs(2).type = 'cylinder';
%obs(2).xy   = [0.60; 0.10];
%obs(2).R    = 0.05;
%obs(2).z    = [0.0; 1.2];

% Initial state and control guesses
xInitial = [x_start; l_start; th_start; phi_start; psi_start; ...
    dx_start; dl_start; dth_start; dphi_start; dpsi_start];
xFinal = [x_end; l_end; th_end; phi_end; psi_end; ...
    dx_end; dl_end; dth_end; dphi_end; dpsi_end];

%% Z Matrix Construction

nx = 10; nu = 3; % number of states and controls

idx_x = reshape(1 : nx*(N+1), nx, N+1);
idx_u = reshape(nx*(N+1) + (1 : nu*(N+1)), nu, N+1);
idx_tf = nx*(N+1) + nu*(N+1) + 1; % add as final state in Z

% Make lb and ub

x_lb = [x_min; l_min; -inf; -inf; -inf; dx_min; dl_min; dth_min; -inf; -inf];
x_ub = [x_max; l_max; inf; inf; inf; dx_max; dl_max; dth_max; inf; inf];

u_lb = [ddx_min; ddl_min; ddth_min];
u_ub = [ddx_max; ddl_max; ddth_max];

lb = zeros(idx_tf, 1);
ub = zeros(size(lb));

for k = 1:(N+1)
    lb(idx_x(:,k)) = x_lb;
    ub(idx_x(:,k)) = x_ub;
    lb(idx_u(:,k)) = u_lb;
    ub(idx_u(:,k)) = u_ub;
end

lb(idx_tf) = tf_min;
ub(idx_tf) = tf_max;

% Pin state start (k=1) and end state (k=N+1)
lb(idx_x(:,1)) = xInitial;
ub(idx_x(:,1)) = xInitial;
lb(idx_x(:,N+1)) = xFinal;
ub(idx_x(:,N+1)) = xFinal;


%% Build initial guess Z0

% The guess passes through start -> waypoints -> finish, in order, spread evenly over
% the nodes. Waypoints are ground-frame payload XY points (payload hanging straight down);
% only trolley radius and slew change, every other state is copied from the start.
wp_XY = [0.76*cos(th_start), 0.76;      % 1: straight out from the start, so it doesn't clip the wall's corner
         0.76*sin(th_start), 0.05];     % 2: just past the wall's outer end (0.65 + R_pay + margin)

points = xInitial;
for j = 1:size(wp_XY, 2)
    x_wp    = xInitial;
    x_wp(1) = norm(wp_XY(:, j));                    % trolley radius above the point
    x_wp(3) = atan2(wp_XY(2, j), wp_XY(1, j));      % slew angle of the point
    points  = [points, x_wp]; %#ok<AGROW>
end
points = [points, xFinal];

Z0 = generateGuess(points, idx_x, numel(lb));

Z0(idx_tf) = tf_guess;
% Linear spread from initial to final

%% fmincon call

options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 1000);

Z_opt = fmincon(@(Z) cost_fun(Z, idx_x, idx_u, idx_tf, N, u_ub, w_t), Z0, ...
    [], [], [], [], lb, ub, ...
    @(Z) nonlcon_fun(Z, idx_x, idx_u, idx_tf, N, nx, obs, H_jib, R_pay), options);

% Unpack optimal trajectories into clean matrices
tf_opt = Z_opt(idx_tf);
t_grid = linspace(0, tf_opt, N+1);
x_opt  = Z_opt(idx_x); % 10 x (N+1)
u_opt  = Z_opt(idx_u); % 3 x (N+1)

%% Save the trajectory for sim_trajectory.m
traj.t       = t_grid;        % 1 x (N+1) node times (s)
traj.x       = x_opt;         % 10 x (N+1) states at the nodes
traj.u       = u_opt;         % 3 x (N+1) inputs at the nodes [ddx; ddl; ddth]
traj.tf      = tf_opt;
traj.w_t     = w_t;           % time weight in the cost (shown on the sim plot)
traj.N       = N;
traj.x_lb    = x_lb;   traj.x_ub = x_ub;      % limits, so the sim can check them
traj.u_lb    = u_lb;   traj.u_ub = u_ub;
traj.obs     = obs;           % obstacle struct array (see Obstacles section)
traj.H_jib   = H_jib;         % jib height used for the obstacle constraints
traj.R_pay   = R_pay;         % payload radius used for the obstacle constraints
traj.created = char(datetime('now'));
out_file = fullfile(fileparts(mfilename('fullpath')), 'trajectory.mat');
save(out_file, '-struct', 'traj');
fprintf('Saved trajectory to %s\n', out_file);

function J = cost_fun(Z, idx_x, idx_u, idx_tf, N, u_ub, w_t)
    x = Z(idx_x);
    u = Z(idx_u);
    tf = Z(idx_tf);
    dt = tf / N;

    R = diag(1 ./ (u_ub.^2)); % Bryson's Rule

    effort = 0.5 * sum((R * u) .* u, 1); % Compute 0.5 * u' * R * u at each time step k


    % penalize residual vibration further in the last 10% of motion
    tail_idx = round(0.9 * (N-1)) : (N-1);
    w_tail = 200.0; %Weight of this punishment 
    tail_pen = w_tail * sum(x(4, tail_idx).^2 + x(5, tail_idx).^2 +x(9, tail_idx).^2 +x(10, tail_idx).^2 );

    % trapezoidal integral sum

    J = sum(0.5 * dt * (effort(1:end-1) + effort(2:end))) + w_t * tf + tail_pen;
end

function [c, ceq] = nonlcon_fun(Z, idx_x, idx_u, idx_tf, N, nx, obs, H_jib, R_pay)
    x = Z(idx_x);
    u = Z(idx_u);
    tf = Z(idx_tf);
    dt = tf / N;

    ceq = zeros(nx*N, 1);

    for k = 1:N
        fk = dynamics(x(:, k), u(:, k));
        fk1 = dynamics(x(:, k+1), u(:, k+1));

        defect = x(:, k+1) - x(:,k) - 0.5*dt*(fk+fk1);
        ceq((k-1)*nx + (1:nx)) = defect;
    end

    num_obs = numel(obs);
    c = zeros(num_obs * (N+1), 1);
    row = 1;

    for k = 1:(N+1)
        % Map states into Payload Position [Xp; Yp; Zp]
        xk = x(1, k); lk = x(2, k); thk = x(3, k);
        phik = x(4,k); psik = x(5,k);

        r_loc = xk + lk *sin(phik)*cos(psik);
        t_loc = lk*sin(psik);
        z_loc = H_jib - lk*cos(phik)*cos(psik);

        Xp = r_loc *cos(thk) - t_loc*sin(thk);
        Yp = r_loc *sin(thk) + t_loc*cos(thk);
        Zp = z_loc;
        p_pay = [Xp; Yp; Zp];

        % Enforce SDF >= 0  ->  -SDF <= 0
        for j = 1:num_obs
            c(row + j - 1) = -obstacle_sdf(p_pay, obs(j), R_pay);
        end
        row = row + num_obs;
    end
end

function Z0 = generateGuess(points, idx_x, nZ)
% Piecewise-linear state guess through the columns of points (nx x P), spread evenly
% over the N+1 nodes. Inputs stay 0; t_f is set by the caller.
    Z0 = zeros(nZ, 1);
    P  = size(points, 2);
    s  = linspace(1, P, size(idx_x, 2));               % node k -> position along the points
    Z0(idx_x) = interp1(1:P, points.', s, 'linear').'; % nx x (N+1)
end