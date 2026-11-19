%% run_sim.m - Plot and animate a saved maneuver (no solving)
% Shows any solution in results/: optimal control (solve_oc.m) or input shaping (solve_is.m).
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'model', 'generated'));

solution = 'oc_min_time.mat';      % or 'oc_fixed.mat', 'is_classic.mat', 'is_cancel.mat', 'is_robust.mat'
fps      = 50;     % animation frame rate (also the plot sampling rate)
t_coast  = 5;      % seconds of U = 0 after the maneuver, to watch for residual swing
animate  = true;   % false: plots only

%% 1. Load
in_file = crane.file(solution);
if ~isfile(in_file), error('No %s yet: run solve_oc.m or solve_is.m first.', in_file); end
S   = load(in_file);
cfg = S.cfg;  t_f = S.t_f;
is_optimal = strcmp(S.kind, 'optimal');
fprintf('Loaded %s\n  made by %s on %s | t_f = %.2f s | loops = %+d\n', in_file, S.source, S.created, t_f, cfg.kdir);
if is_optimal, fprintf('  optimal control, time price W_t = %.3g\n', cfg.W_t);
else,          fprintf('  input-shaped command, simulated on the full nonlinear model\n'); end

%% 2. Evaluate the maneuver
t_man = linspace(0, t_f, round(t_f*fps));
T = crane.trajectory(S, t_man);
M = crane.metrics(T, cfg);

%% 3. Coast: freeze the actuators (U = 0) and integrate the full nonlinear model forward
stopped = struct('U', @(t) zeros(3, numel(t)), 'breaks', [], 'pwc', false);
Tc  = crane.simulate(stopped, t_f + (0 : 1/fps : t_coast), cfg, T.X(:, end));
Xc  = Tc.X(:, 2:end);  t_c = Tc.t(2:end);   % drop the repeated t_f point
Gc  = structfun(@(v) v(:, 2:end), Tc.G, 'UniformOutput', false);

% Maneuver + coast, in the variable names the plots use
k_end  = numel(t_man);                      % index of the last maneuver sample (t = t_f)
t_sim  = [t_man, t_c];
Nt     = numel(t_sim);
X_traj = [T.X, Xc];
U_traj = [T.U, zeros(3, numel(t_c))];
H_t    = [T.H, nan(1, numel(t_c))];         % H only means something under the optimal control
Qp     = [T.G.Qp, Gc.Qp];
r_rad  = [T.G.r_rad, Gc.r_rad];
d_pin  = [T.G.d_pin, Gc.d_pin];
d_cable = [T.G.d_cable, Gc.d_cable];

x_t   = X_traj(1, :);  phi_t = X_traj(4, :);  psi_t = X_traj(5, :);
dx_t  = X_traj(6, :);  dl_t  = X_traj(7, :);  dth_t = X_traj(8, :);
eta_t = X_traj(11, :);

xmax = cfg.xmax;  R_p = cfg.R_p;  r_p = cfg.r_p;  c_clr = cfg.c_clr;  z_top = cfg.z_top;
H_tower = cfg.lmax + 0.15;                  % jib height above the floor (drawing only)
lift  = [0; 0; H_tower];
Jib_I = [T.G.Jib, Gc.Jib] + lift;
P_I   = [T.G.P,   Gc.P]   + lift;
Q_I   = Qp + lift;

%% 4. Checks
fprintf('\n--- Solution checks ---\n');
fprintf('final time t_f          : %.2f s\n', t_f);
if is_optimal, fprintf('H spread (max-min)      : %.2e   (should be ~0)\n', M.H_spread); end
fprintf('max x / xmax            : %.3f   (must be <= 1)\n', M.x_ratio);
fprintf('max |v| / vmax (x,l,th) : %.3f %.3f %.3f\n', M.v_ratio);
fprintf('max |U| / Umax (x,l,th) : %.3f %.3f %.3f\n', M.u_ratio);
fprintf('min payload-pin dist    : %.1f mm (clearance %.1f mm)\n', 1000*M.d_pin_min, 1000*(r_p + c_clr));
fprintf('min cable-pin dist      : %.1f mm (pin radius %.1f mm)\n', 1000*M.d_cable_min, 1000*r_p);
fprintf('max payload reach - xmax: %+.1f mm  (> 0 means it swung past the trolley limit)\n', 1000*M.reach_over);
fprintf('windings eta(tf)/2pi    : %.4f\n', M.windings);
fprintf('swing at t_f phi, psi   : %.3f, %.3f deg\n', M.end_swing);
fprintf('leftover swing at t_f   : %.1f mm (radial %.1f, tangential %.1f)\n', 1000*M.residual, 1000*M.residual_parts);
if M.admissible, fprintf('ADMISSIBLE: all limits respected, pin never touched\n');
else,            fprintf('NOT ADMISSIBLE: %s\n', strjoin(M.why, ', ')); end

coast = (k_end+1):Nt;
fprintf('\n--- Coast (%g s, U = 0) ---\n', t_coast);
fprintf('max |phi|, |psi|        : %.2e, %.2e deg\n', rad2deg(max(abs(phi_t(coast)))), rad2deg(max(abs(psi_t(coast)))));
fprintf('max payload drift       : %.3f mm from its position at t_f\n', 1000*max(vecnorm(Qp(:, coast) - Qp(:, k_end))));

%% 5. Telemetry (mm and deg)
figure('Name', sprintf('Telemetry: %s', S.source), 'Color', 'w', 'Position', [50, 50, 1150, 750]);

subplot(2,2,1); hold on; grid on;
plot(t_sim, x_t*1000,   'b-', 'LineWidth', 2, 'DisplayName', 'Trolley x(t)');
plot(t_sim, r_rad*1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Payload reach');
yline(xmax*1000, 'k--', 'LineWidth', 1.5, 'DisplayName', 'x_{max} limit');
yline((R_p + r_p + c_clr)*1000, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Reach needed to clear pin');
xlabel('Time (s)'); ylabel('Radius (mm)');
title('Payload Reach vs. Trolley Limit'); legend('Location', 'best');

subplot(2,2,2); hold on; grid on;
plot(t_sim, rad2deg(phi_t), 'r-', 'LineWidth', 2, 'DisplayName', '\phi (radial swing)');
plot(t_sim, rad2deg(psi_t), 'm-', 'LineWidth', 2, 'DisplayName', '\psi (tangential swing)');
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Angle (deg)');
title(sprintf('Swing Angles (at t_f: \\phi=%.2f^\\circ, \\psi=%.2f^\\circ)', M.end_swing));
legend('Location', 'best');

subplot(2,2,3); hold on; grid on;
plot(t_sim, dx_t*1000,      'b-', 'LineWidth', 1.8, 'DisplayName', 'dx (max 72.6 mm/s)');
plot(t_sim, dl_t*1000,      'g-', 'LineWidth', 1.8, 'DisplayName', 'dl (max 30.0 mm/s)');
plot(t_sim, rad2deg(dth_t), 'c-', 'LineWidth', 1.8, 'DisplayName', 'd\theta (max 32.4 deg/s)');
yline([72.6, -72.6], 'b:', 'HandleVisibility', 'off');
yline([30, -30],     'g:', 'HandleVisibility', 'off');
yline([32.4, -32.4], 'c:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Velocity (mm/s or deg/s)');
title('Actuator Velocities vs. Hardware Limits'); legend('Location', 'best');

subplot(2,2,4); hold on; grid on;
plot(t_sim, U_traj(1,:)*1000,     'b-', 'LineWidth', 1.8, 'DisplayName', 'ddx (max 72.6 mm/s^2)');
plot(t_sim, U_traj(2,:)*1000,     'g-', 'LineWidth', 1.8, 'DisplayName', 'ddl (max 30.0 mm/s^2)');
plot(t_sim, rad2deg(U_traj(3,:)), 'c-', 'LineWidth', 1.8, 'DisplayName', 'dd\theta (max 32.4 deg/s^2)');
yline([72.6, -72.6], 'b:', 'HandleVisibility', 'off');
yline([30, -30],     'g:', 'HandleVisibility', 'off');
yline([32.4, -32.4], 'c:', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Accel (mm/s^2 or deg/s^2)');
title('Commanded Accelerations'); legend('Location', 'best');

mark_coast(gcf, t_f);

%% 6. Obstacle & optimality checks
figure('Name', sprintf('Pin Loop Checks: %s', S.source), 'Color', 'w', 'Position', [100, 80, 1150, 750]);
ang = linspace(0, 2*pi, 200);

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

subplot(2,2,2); hold on; grid on;
plot(t_sim, d_pin*1000,   'r-', 'LineWidth', 1.8, 'DisplayName', 'Payload to pin center');
plot(t_sim, d_cable*1000, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Cable (below pin top) to pin center');
yline((r_p + c_clr)*1000, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Payload clearance');
yline(r_p*1000, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Pin radius');
xlabel('Time (s)'); ylabel('Distance (mm)'); title('Pin Clearance'); legend('Location', 'best');

subplot(2,2,3); hold on; grid on;
if is_optimal
    plot(t_sim, H_t, 'LineWidth', 1.8);
    xlabel('Time (s)'); ylabel('H');
    title(sprintf('Hamiltonian (spread %.1e, should be flat)', M.H_spread));
else
    % No Hamiltonian for a shaped command: show the swing it builds up and then cancels
    offset = vecnorm(Qp(1:2,:) - [T.G.P(1:2,:), Gc.P(1:2,:)]);   % payload offset from the trolley
    plot(t_sim, offset*1000, 'LineWidth', 1.8);
    xlabel('Time (s)'); ylabel('Offset (mm)');
    title(sprintf('Payload Swing Offset (left at t_f: %.1f mm)', 1000*M.residual));
end

subplot(2,2,4); hold on; grid on;
plot(t_sim, eta_t/(2*pi), 'LineWidth', 1.8);
yline(cfg.kdir, 'k--', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\eta / 2\pi');
title('Windings Around the Pin');

mark_coast(gcf, t_f);

%% 7. 3D animation (mm)
if ~animate, return; end
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

[Xcyl, Ycyl, Zcyl] = cylinder(r_p, 40);   % unit-height tube, radius r_p, centered at the origin
Zcyl = Zcyl * (H_tower + z_top);          % stretch: floor up to the pin top
surf(ax_3d, (Xcyl + R_p)*1000, Ycyl*1000, Zcyl*1000, ...
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
function mark_coast(fig, t_f)
% Draw a vertical line at t_f on every time-axis plot in the figure.
    for ax = findobj(fig, 'Type', 'axes')'
        if strcmp(ax.XLabel.String, 'Time (s)')
            xline(ax, t_f, 'k-.', 'coast', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end
end
