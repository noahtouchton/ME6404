%% Tower Crane 3D Optimal Control Simulator (With Exact Lab Crane Limits)
clear; clc; close all;

%% 1. Define Physical Parameters & Exact Hardware Limits (in Meters & Radians)
g       = 9.81;                % Gravity (m/s^2)

% Workspace limits (Adjust xmax and lmax to match your physical crane dimensions)
xmax    = 0.925;               % Max trolley position (450 mm = 0.45 m)
lmax    = 1.450;               % Max hoist cable length (600 mm = 0.60 m)

% Exact Velocity Limits:
dxmax   = 72.6e-3;             % 72.6 mm/s  -> 0.0726 m/s
dlmax   = 30.0e-3;             % 30.0 mm/s  -> 0.0300 m/s
dthmax  = deg2rad(32.4);       % 32.4 deg/s -> 0.5655 rad/s

% Exact Acceleration Limits (Same numerical magnitude as velocities!):
ddxmax  = 72.6e-3;             % 72.6 mm/s^2  -> 0.0726 m/s^2
ddlmax  = 30.0e-3;             % 30.0 mm/s^2  -> 0.0300 m/s^2
ddthmax = deg2rad(32.4);       % 32.4 deg/s^2 -> 0.5655 rad/s^2

% Pin

r_p = 0.03;  % Pin radius (m)
c_clr = 0.02; % Clearance: payload size + safety margin (m)
z_top = -0.60; % Pin top relative to jib height (m)

W_obs = 5; % Obstacle wall strength
k_x = 1; % Trolley-limit wall strength
W_t = 0; % Time price


% Optimal Control Weights
r1      = 0.2;                 % Trolley effort weight
r2      = 0.5;                 % Hoist effort weight
r3      = 0.2;                 % Slew effort weight

% Start = end state

x0 = 0.70; l0 = 1.00;

X0 = [x0; l0; 0; 0; 0; 0; 0; 0; 0; 0];
kdir = 1; % +1 is one counter clockwise loop around the pin. -1 for clockwise
t_f = 30;   % warm-up horizon. Longer is NOT easier: 25 s and 40 s failed for 2 loops (costates get stiff)

R_p_start = 0.7;
R_p_final = 0.9;

make_params = @(R_p, k_x, W_t) [g; r1; r2; r3; xmax; lmax; dxmax; dlmax; dthmax; ...
    ddxmax; ddlmax; ddthmax; R_p; 0; r_p; c_clr; W_obs; k_x; W_t];


% Initial Guess: a loop around the pin with the payload hanging straight
% down

opts = bvpset('RelTol', 1e-3, 'AbsTol', 1e-4, 'Nmax', 20000);
T_sw = 2*pi*sqrt(l0/g);
t_mesh = linspace(0, t_f, round(20*t_f/T_sw));
pin0 = [R_p_start; 0];
sol = bvpinit(t_mesh, @(t) loop_guess(t, t_f, X0, pin0, kdir));
% Costates exactly 0 make every lambda*df/dX term vanish -> singular Newton matrix.
% A tiny random costate guess breaks that degeneracy (rng fixed so runs repeat).
rng(1);
sol.y(12:22, :) = 1e-3 * randn(11, numel(sol.x));

% First solve at the easiest pin location
[sol, ok, info] = try_solve(sol, make_params(R_p_start, k_x, W_t), X0, kdir, opts);
fprintf('R_p = %.4f  %s\n', R_p_start, info);
if ~ok, error('First solve failed. Try larger t_f, smaller W_obs, or smaller R_p_start.'); end

% Move the pin outward until swinging is required
R_p  = R_p_start;
step = 0.01;
while R_p < R_p_final - 1e-9
    R_try = min(R_p + step, R_p_final);
    [sol_try, ok, info] = try_solve(sol, make_params(R_try, k_x, W_t), X0, kdir, opts);
    fprintf('R_p = %.4f  %s\n', R_try, info);
    if ok
        sol  = sol_try;                 % accept: it becomes the next guess
        R_p  = R_try;
        step = min(1.5*step, 0.02);     % things are going well: take bigger steps
    else
        step = step/2;                  % failed: retry from the last good sol with a smaller step
        if step < 1e-4, error('Stalled at R_p = %.4f', R_p); end
    end
end

% Harden the trolley wall until the trolley stays inside xmax
while max(deval(sol, linspace(0, t_f, 2000), 1)) > xmax
    k_try = 3*k_x;
    [sol_try, ok, info] = try_solve(sol, make_params(R_p, k_try, W_t), X0, kdir, opts);
    fprintf('k_x = %.0f  %s\n', k_try, info);
    if ~ok, error('Could not harden the x wall.'); end
    sol = sol_try;  k_x = k_try;
end
params = make_params(R_p, k_x, W_t);

% Free final time: trade effort for speed
Z_end = deval(sol, t_f);
H_fix = hamiltonian_value(Z_end, params);            % H of the fixed-tf solution (should be < 0)
fprintf('Fixed-tf H = %.4g  ->  starting W_t = %.4g\n', H_fix, -H_fix);

solT       = [];                                     % stays empty if no free-tf solve succeeds
solT_guess = bvpinit(sol.x / t_f, @(tau) deval(sol, tau*t_f), t_f);   % same solution on tau in [0,1]
for W_try = -H_fix * [1, 1.5, 2.25, 3.4, 5]
    p_try = make_params(R_p, k_x, W_try);
    odeT  = @(tau, Z, T) T * hamiltonian_ode(T*tau, Z, p_try);
    bcT   = @(ya, yb, T) [bc_loop(ya, yb, X0, kdir); hamiltonian_value(yb, p_try)];
    lastwarn('');
    try
        solT_try = bvp4c(odeT, bcT, solT_guess, opts);
    catch
        break
    end
    if ~isempty(lastwarn) || solT_try.stats.maxres > 1e-3, break; end
    solT = solT_try;  solT_guess = solT_try;  W_t = W_try;  params = p_try;
    fprintf('W_t = %.3g  ->  t_f = %.2f s\n', W_t, solT.parameters);
end

if isempty(solT)
    fprintf('Free-tf solve did not converge; keeping fixed t_f = %.2f s\n', t_f);
    sol_eval = @(t) deval(sol, t);
else
    t_f = solT.parameters;
    sol_eval = @(t) deval(solT, t / t_f);            % evaluate in real seconds
end

%% 4. Evaluate Trajectory & Kinematics in 3D Inertial Frame
fps    = 50;
t_sim  = linspace(0, t_f, round(t_f * fps));
Nt     = numel(t_sim);
XP     = sol_eval(t_sim);                 % 22 x Nt
X_traj = XP(1:11, :);

U_traj = zeros(3, Nt);                    % optimal accelerations
H_t    = zeros(1, Nt);                    % Hamiltonian (should be flat)
Qp     = zeros(3, Nt);                    % payload [X; Y; Z] in ground frame, Z relative to jib
for k = 1:Nt
    U_traj(:, k) = optimal_inputs(XP(:, k), params);
    H_t(k)       = hamiltonian_value(XP(:, k), params);
    Qp(:, k)     = payload_position(X_traj(:, k), params);
end

% Coast: freeze the actuators (U = 0) for t_coast seconds after the maneuver
% and integrate the real nonlinear EOM forward to watch for residual swing.
% hamiltonian_ode with all costates = 0 gives U* = 0, so its first 11 rows are f(X, U=0).
t_coast = 5;
t_c     = t_f + (0 : 1/fps : t_coast);
coast_f = @(t, X) first11(hamiltonian_ode(t, [X; zeros(11,1)], params));
[~, Xc] = ode45(coast_f, t_c, X_traj(:, end), odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
Xc  = Xc(2:end, :).';                     % ode45 returns rows; drop the repeated t_f point
t_c = t_c(2:end);
Qc  = zeros(3, numel(t_c));
for k = 1:numel(t_c)
    Qc(:, k) = payload_position(Xc(:, k), params);
end

k_end  = Nt;                              % index of the last maneuver sample (t = t_f)
t_sim  = [t_sim, t_c];                    % from here on, every array covers maneuver + coast
X_traj = [X_traj, Xc];
U_traj = [U_traj, zeros(3, numel(t_c))];
H_t    = [H_t, nan(1, numel(t_c))];       % H only means something under the optimal control
Qp     = [Qp, Qc];
Nt     = numel(t_sim);

x_t   = X_traj(1, :);  l_t   = X_traj(2, :);  th_t  = X_traj(3, :);
phi_t = X_traj(4, :);  psi_t = X_traj(5, :);
dx_t  = X_traj(6, :);  dl_t  = X_traj(7, :);  dth_t = X_traj(8, :);
eta_t = X_traj(11, :);

H_tower = lmax + 0.15;                    % jib height above the floor (drawing only)
Jib_I = [xmax*cos(th_t); xmax*sin(th_t); H_tower*ones(1, Nt)];   % jib tip
P_I   = [x_t.*cos(th_t);  x_t.*sin(th_t);  H_tower*ones(1, Nt)];  % trolley
Q_I   = Qp + [0; 0; H_tower];                                     % payload

r_rad = hypot(Qp(1,:), Qp(2,:));          % payload radius from the tower
d_pin = hypot(Qp(1,:) - R_p, Qp(2,:));    % payload distance to the pin center

% Cable vs pin: the part of the cable below the pin top must stay outside the pin
d_cable = inf(1, Nt);
for k = 1:Nt
    if Qp(3,k) < z_top                                   % payload hangs below the pin top
        s_top = z_top / Qp(3,k);                         % fraction of the cable above pin-top height
        A     = P_I(1:2,k) + s_top*(Qp(1:2,k) - P_I(1:2,k));   % cable point at pin-top height
        d_cable(k) = seg_point_dist(A, Qp(1:2,k), [R_p; 0]);
    end
end

fprintf('\n--- Solution checks ---\n');
fprintf('final time t_f          : %.2f s\n', t_f);
fprintf('H spread (max-min)      : %.2e   (should be ~0)\n', max(H_t) - min(H_t));
fprintf('max x / xmax            : %.3f   (must be <= 1)\n', max(x_t)/xmax);
fprintf('max |dx| / dxmax        : %.3f\n', max(abs(dx_t))/dxmax);
fprintf('max |dl| / dlmax        : %.3f\n', max(abs(dl_t))/dlmax);
fprintf('max |dth| / dthmax      : %.3f\n', max(abs(dth_t))/dthmax);
fprintf('max |U| / Umax          : %.3f %.3f %.3f\n', max(abs(U_traj), [], 2) ./ [ddxmax; ddlmax; ddthmax]);
fprintf('min payload-pin dist    : %.1f mm (clearance %.1f mm)\n', 1000*min(d_pin), 1000*(r_p + c_clr));
fprintf('min cable-pin dist      : %.1f mm (pin radius %.1f mm)\n', 1000*min(d_cable), 1000*r_p);
fprintf('max payload reach - xmax: %+.1f mm  (> 0 means it swung past the trolley limit)\n', 1000*(max(r_rad) - xmax));
fprintf('windings eta(tf)/2pi    : %.4f\n', eta_t(end)/(2*pi));
fprintf('swing at t_f phi, psi   : %.3f, %.3f deg\n', rad2deg(phi_t(k_end)), rad2deg(psi_t(k_end)));

coast = (k_end+1):Nt;
fprintf('\n--- Coast (%g s, U = 0) ---\n', t_coast);
fprintf('max |phi|, |psi|        : %.2e, %.2e deg\n', rad2deg(max(abs(phi_t(coast)))), rad2deg(max(abs(psi_t(coast)))));
fprintf('max payload drift       : %.3f mm from its position at t_f\n', 1000*max(vecnorm(Qp(:, coast) - Qp(:, k_end))));
fprintf('max |dx|, |dl|, |dth|   : %.2e mm/s, %.2e mm/s, %.2e deg/s\n', ...
    1000*max(abs(dx_t(coast))), 1000*max(abs(dl_t(coast))), rad2deg(max(abs(dth_t(coast)))));

%% 5. Plot Telemetry (mm and deg)
figure('Name', 'Hardware-Constrained Optimal Control', 'Color', 'w', 'Position', [50, 50, 1150, 750]);

% Subplot 1: Radial reach (mm)
subplot(2,2,1); hold on; grid on;
plot(t_sim, x_t*1000,   'b-', 'LineWidth', 2, 'DisplayName', 'Trolley x(t)');
plot(t_sim, r_rad*1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Payload reach');
yline(xmax*1000, 'k--', 'LineWidth', 1.5, 'DisplayName', 'x_{max} limit');
yline((R_p + r_p + c_clr)*1000, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Reach needed to clear pin');
xlabel('Time (s)'); ylabel('Radius (mm)');
title('Payload Reach vs. Trolley Limit'); legend('Location', 'best');

% Subplot 2: Swing angles (deg)
subplot(2,2,2); hold on; grid on;
plot(t_sim, rad2deg(phi_t), 'r-', 'LineWidth', 2, 'DisplayName', '\phi (radial swing)');
plot(t_sim, rad2deg(psi_t), 'm-', 'LineWidth', 2, 'DisplayName', '\psi (tangential swing)');
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Angle (deg)');
title(sprintf('Swing Angles (final: \\phi=%.2f^\\circ, \\psi=%.2f^\\circ)', rad2deg(phi_t(end)), rad2deg(psi_t(end))));
legend('Location', 'best');

% Subplot 3: Actuator velocities vs limits
subplot(2,2,3); hold on; grid on;
plot(t_sim, dx_t*1000,      'b-', 'LineWidth', 1.8, 'DisplayName', 'dx (max 72.6 mm/s)');
plot(t_sim, dl_t*1000,      'g-', 'LineWidth', 1.8, 'DisplayName', 'dl (max 30.0 mm/s)');
plot(t_sim, rad2deg(dth_t), 'c-', 'LineWidth', 1.8, 'DisplayName', 'd\theta (max 32.4 deg/s)');
yline([72.6, -72.6], 'b:', 'HandleVisibility', 'off');
yline([30, -30],     'g:', 'HandleVisibility', 'off');
yline([32.4, -32.4], 'c:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Velocity (mm/s or deg/s)');
title('Actuator Velocities vs. Hardware Limits'); legend('Location', 'best');

% Subplot 4: Optimal accelerations vs limits
subplot(2,2,4); hold on; grid on;
plot(t_sim, U_traj(1,:)*1000,      'b-', 'LineWidth', 1.8, 'DisplayName', 'ddx (max 72.6 mm/s^2)');
plot(t_sim, U_traj(2,:)*1000,      'g-', 'LineWidth', 1.8, 'DisplayName', 'ddl (max 30.0 mm/s^2)');
plot(t_sim, rad2deg(U_traj(3,:)),  'c-', 'LineWidth', 1.8, 'DisplayName', 'dd\theta (max 32.4 deg/s^2)');
yline([72.6, -72.6], 'b:', 'HandleVisibility', 'off');
yline([30, -30],     'g:', 'HandleVisibility', 'off');
yline([32.4, -32.4], 'c:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Accel (mm/s^2 or deg/s^2)');
title('Optimal Accelerations'); legend('Location', 'best');

mark_coast(gcf, t_f);

%% 5b. Obstacle & optimality checks
figure('Name', 'Pin Loop & Optimality Checks', 'Color', 'w', 'Position', [100, 80, 1150, 750]);
ang = linspace(0, 2*pi, 200);

% Top-down view: the loop around the pin
subplot(2,2,1); hold on; grid on; axis equal;
fill((R_p + r_p*cos(ang))*1000, r_p*sin(ang)*1000, [0.2 0.6 0.2], 'EdgeColor', 'none', 'DisplayName', 'Pin');
plot((R_p + (r_p + c_clr)*cos(ang))*1000, (r_p + c_clr)*sin(ang)*1000, 'g--', 'DisplayName', 'Clearance');
plot(xmax*1000*cos(ang), xmax*1000*sin(ang), 'k--', 'DisplayName', 'x_{max}');
plot(P_I(1,:)*1000, P_I(2,:)*1000, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trolley');
plot(Qp(1,:)*1000,  Qp(2,:)*1000,  'r-', 'LineWidth', 1.5, 'DisplayName', 'Payload');
out = r_rad > xmax;                                   % samples where only swing explains the reach
plot(Qp(1,out)*1000, Qp(2,out)*1000, 'm.', 'MarkerSize', 10, 'DisplayName', 'Beyond x_{max}');
plot(Qp(1,1)*1000, Qp(2,1)*1000, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Start / end');
allX = [Qp(1,:), P_I(1,:), R_p + r_p + c_clr];  allY = [Qp(2,:), P_I(2,:)];
xlim([min(allX) - 0.05, max(allX) + 0.05]*1000);
ylim([min(allY) - 0.05, max(allY) + 0.05]*1000);
xlabel('X (mm)'); ylabel('Y (mm)'); title('Top-Down View'); legend('Location', 'best');

% Distance to the pin
subplot(2,2,2); hold on; grid on;
plot(t_sim, d_pin*1000, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Payload to pin center');
plot(t_sim, d_cable*1000, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Cable (below pin top) to pin center');
yline((r_p + c_clr)*1000, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Payload clearance');
yline(r_p*1000, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Pin radius');
xlabel('Time (s)'); ylabel('Distance (mm)'); title('Pin Clearance'); legend('Location', 'best');

% Hamiltonian
subplot(2,2,3); hold on; grid on;
plot(t_sim, H_t, 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('H');
title(sprintf('Hamiltonian (spread %.1e, should be flat)', max(H_t) - min(H_t)));

% Winding around the pin
subplot(2,2,4); hold on; grid on;
plot(t_sim, eta_t/(2*pi), 'LineWidth', 1.8);
yline(kdir, 'k--', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\eta / 2\pi');
title('Windings Around the Pin');

mark_coast(gcf, t_f);

%% 6. Interactive 3D Animation (mm)
fig_anim = figure('Name', '3D Tower Crane Simulation', 'Color', 'w', 'Position', [200, 100, 900, 700]);
ax_3d = axes('Parent', fig_anim);
hold(ax_3d, 'on'); grid(ax_3d, 'on'); axis(ax_3d, 'equal');
view(ax_3d, 45, 25);
xlim(ax_3d, [-1.3*xmax, 1.3*xmax]*1000);
ylim(ax_3d, [-1.3*xmax, 1.3*xmax]*1000);
zlim(ax_3d, [0, H_tower + 0.1]*1000);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

plot3(ax_3d, xmax*1000*cos(ang), xmax*1000*sin(ang), zeros(size(ang)), 'k--', 'LineWidth', 1.2, 'DisplayName', 'x_{max} limit');
plot3(ax_3d, [0, 0], [0, 0], [0, H_tower*1000], 'Color', [0.3 0.3 0.3], 'LineWidth', 6, 'DisplayName', 'Tower');

[Xc, Yc, Zc] = cylinder(r_p, 40);         % unit-height tube, radius r_p, centered at the origin
Zc = Zc * (H_tower + z_top);              % stretch: floor up to the pin top
surf(ax_3d, (Xc + R_p)*1000, Yc*1000, Zc*1000, ...
     'FaceColor', [0.2 0.6 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'Pin');

h_jib     = plot3(ax_3d, [0, Jib_I(1,1)]*1000, [0, Jib_I(2,1)]*1000, [H_tower, Jib_I(3,1)]*1000, 'Color', [0.85 0.6 0.1], 'LineWidth', 4, 'DisplayName', 'Jib');
h_trolley = plot3(ax_3d, P_I(1,1)*1000, P_I(2,1)*1000, P_I(3,1)*1000, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Trolley');
h_cable   = plot3(ax_3d, [P_I(1,1), Q_I(1,1)]*1000, [P_I(2,1), Q_I(2,1)]*1000, [P_I(3,1), Q_I(3,1)]*1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Cable');
h_payload = plot3(ax_3d, Q_I(1,1)*1000, Q_I(2,1)*1000, Q_I(3,1)*1000, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Payload');
h_trail   = plot3(ax_3d, Q_I(1,1)*1000, Q_I(2,1)*1000, Q_I(3,1)*1000, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Payload path');

legend(ax_3d, 'Location', 'northeast');
title_handle = title(ax_3d, '');

for k = 1:Nt
    if ~isvalid(fig_anim), break; end
    set(h_jib,     'XData', [0, Jib_I(1,k)]*1000, 'YData', [0, Jib_I(2,k)]*1000, 'ZData', [H_tower, Jib_I(3,k)]*1000);
    set(h_trolley, 'XData', P_I(1,k)*1000,        'YData', P_I(2,k)*1000,        'ZData', P_I(3,k)*1000);
    set(h_cable,   'XData', [P_I(1,k), Q_I(1,k)]*1000, 'YData', [P_I(2,k), Q_I(2,k)]*1000, 'ZData', [P_I(3,k), Q_I(3,k)]*1000);
    set(h_payload, 'XData', Q_I(1,k)*1000,        'YData', Q_I(2,k)*1000,        'ZData', Q_I(3,k)*1000);
    set(h_trail,   'XData', Q_I(1,1:k)*1000,      'YData', Q_I(2,1:k)*1000,      'ZData', Q_I(3,1:k)*1000);
    if k > k_end, phase = 'COAST (U = 0)'; else, phase = 'maneuver'; end
    set(title_handle, 'String', sprintf('%s  t = %.2f s (t_f = %.2f) | reach %.1f mm | pin dist %.1f mm | windings %.2f', ...
        phase, t_sim(k), t_f, r_rad(k)*1000, d_pin(k)*1000, eta_t(k)/(2*pi)));
    drawnow;
    pause(1/fps);
end

%% Local functions
function Y = loop_guess(t, t_f, X0, pin, kdir)
% Initial guess: payload travels a circle around the pin, cable straight down.
    tau = t / t_f;
    s   = tau - sin(2*pi*tau)/(2*pi);          % 0 -> 1, zero slope at both ends (starts/stops at rest)
    sd  = (1 - cos(2*pi*tau)) / t_f;           % ds/dt
    rho = norm([X0(1); 0] - pin);              % loop radius = start's distance from the pin
    a   = 2*pi*abs(kdir)*s;                    % angle around the pin: |kdir| full turns
    Q   = pin + rho*[-cos(a); -sign(kdir)*sin(a)];            % at a = 0 this is exactly the start point
    Qd  = rho*2*pi*abs(kdir)*sd*[sin(a); -sign(kdir)*cos(a)]; % its time derivative

    x   = norm(Q);                             % trolley radius (directly above the payload)
    th  = atan2(Q(2), Q(1));                   % jib angle
    dx  = dot(Q, Qd) / x;                      % radial rate
    dth = (Q(1)*Qd(2) - Q(2)*Qd(1)) / x^2;     % angular rate about the tower
    eta = kdir*2*pi*s;                         % winding angle so far

    Y = [x; X0(2); th; 0; 0; dx; 0; dth; 0; 0; eta; zeros(11,1)];   % costate guess = 0
end

function res = bc_loop(ya, yb, X0, ~)       % 4th input (kdir) no longer used: see below
% 22 boundary conditions: start at X0, eta = 0; end at X0; p_eta = 0.
% eta(tf) = 2*pi*kdir is NOT imposed: eta is a function of the path's endpoints
% and topology, so that BC is degenerate and leaves p_eta undetermined (singular
% Jacobian). The loop is set by the guess + the pin wall; try_solve checks eta(tf).
    res = [ya(1:10) - X0;
           ya(11);
           yb(1:10) - X0;
           yb(22)];
end

function [sol, ok, info] = try_solve(sol_guess, params, X0, kdir, opts)
% One bvp4c solve with every failure check. Returns the guess unchanged on failure.
    odefun = @(t, Z) hamiltonian_ode(t, Z, params);
    bcfun  = @(ya, yb) bc_loop(ya, yb, X0, kdir);
    lastwarn('');                                   % clear old warnings
    try
        sol = bvp4c(odefun, bcfun, sol_guess, opts);
    catch err                                       % e.g. "singular Jacobian"
        sol = sol_guess;  ok = false;  info = err.message;  return
    end
    msg = lastwarn;                                 % non-empty = bvp4c complained
    Zf  = sol.y(:, end);
    end_err = max(abs([Zf(1:10) - X0; Zf(11) - 2*pi*kdir]));
    ok   = isempty(msg) && sol.stats.maxres < 1e-3 && end_err < 1e-3;
    info = sprintf('maxres=%.1e end_err=%.1e mesh=%d %s', ...
                   sol.stats.maxres, end_err, numel(sol.x), msg);
    if ~ok, sol = sol_guess; end
end

function d = first11(F)
% Keep the 11 state derivatives of the 22-row Hamiltonian ODE.
    d = F(1:11);
end

function mark_coast(fig, t_f)
% Draw a vertical line at t_f on every time-axis plot in the figure.
    for ax = findobj(fig, 'Type', 'axes')'
        if strcmp(ax.XLabel.String, 'Time (s)')
            xline(ax, t_f, 'k-.', 'coast', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end
end

function d = seg_point_dist(A, B, C)
% Shortest distance from point C to the line segment A-B (all 2x1).
    AB = B - A;
    s  = max(0, min(1, dot(C - A, AB) / max(dot(AB, AB), eps)));   % closest point's fraction along AB
    d  = norm(A + s*AB - C);
end
