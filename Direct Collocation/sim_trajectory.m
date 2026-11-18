%% sim_trajectory.m - Run the crane sim with the trajectory saved by optimize.m
% Loads trajectory.mat, drives the nonlinear dynamics (dynamics.m) with the optimal inputs,
% then lets it coast with u = 0 to show any residual swing. Compares the simulated states
% with the optimizer's nodes, checks the limits, plots and animates.
%
% u(t) is interpolated LINEARLY between nodes, because that's what trapezoidal collocation
% assumes: x(k+1) = x(k) + dt/2·(f(x_k, u_k) + f(x_k+1, u_k+1)).
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'trajectory.mat'));

t_coast = 5;      % s of u = 0 after the maneuver
animate = true;   % 3D animation at the end
fps     = 30;

% Obstacles saved by optimize.m (older trajectory files have none)
has_obs = isfield(S, 'obs') && ~isempty(S.obs);
if isfield(S, 'H_jib'), H = S.H_jib; else, H = max(S.x_ub(2), 1) + 0.2; end   % jib height above the floor
if isfield(S, 'R_pay'), R_pay = S.R_pay; else, R_pay = 0; end

%% 1. Simulate: maneuver, then coast
u_of_t = @(t) interp1(S.t, S.u.', t, 'linear', 0).';     % 3x1, zero outside [0, tf]
opts   = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', S.tf/S.N/4);

[t1, X1] = ode45(@(t, x) dynamics(x, u_of_t(t)), [0 S.tf], S.x(:, 1), opts);
[t2, X2] = ode45(@(t, x) dynamics(x, zeros(3, 1)), [S.tf, S.tf + t_coast], X1(end, :).', opts);

t = [t1; t2(2:end)];   X = [X1; X2(2:end, :)].';           % X is 10 x Nt
U = [u_of_t(t1.'), zeros(3, numel(t2) - 1)];

%% 2. Checks
Xn = interp1(t, X.', S.t).';                                % sim at the optimizer's nodes
err = abs(Xn - S.x);
names = {'x', 'l', 'th', 'phi', 'psi', 'dx', 'dl', 'dth', 'dphi', 'dpsi'};
fprintf('Trajectory from %s: N = %d, t_f = %.2f s\n\n', S.created, S.N, S.tf);
fprintf('Max |sim - optimizer| over the nodes (collocation error):\n');
for i = 1:10, fprintf('  %-5s %.2e\n', names{i}, max(err(i, :))); end

xf_err = X1(end, :).' - S.x(:, end);
fprintf('\nEnd of maneuver, sim - target: x %+.2f mm, l %+.2f mm, th %+.3f deg, phi %+.3f deg, psi %+.3f deg\n', ...
        1000*xf_err(1), 1000*xf_err(2), rad2deg(xf_err(3)), rad2deg(xf_err(4)), rad2deg(xf_err(5)));

coast = t > S.tf;
Q  = payload_xyz(X);
Qf = Q(:, find(~coast, 1, 'last'));
fprintf('Residual swing during the coast: max payload drift %.2f mm, max |phi| %.3f deg, max |psi| %.3f deg\n', ...
        1000*max(vecnorm(Q(1:2, coast) - Qf(1:2))), rad2deg(max(abs(X(4, coast)))), rad2deg(max(abs(X(5, coast)))));

viol = {};
for i = 1:10
    if any(X(i, ~coast) < S.x_lb(i) - 1e-6 | X(i, ~coast) > S.x_ub(i) + 1e-6), viol{end+1} = names{i}; end %#ok<SAGROW>
end
if any(any(U(:, ~coast) < S.u_lb - 1e-9 | U(:, ~coast) > S.u_ub + 1e-9)), viol{end+1} = 'u'; end
if isempty(viol), fprintf('Limits: all respected\n'); else, fprintf('Limits VIOLATED: %s\n', strjoin(viol, ', ')); end

if has_obs
    Qg = Q + [0; 0; H];                                     % payload in the ground frame (Z up from floor)
    fprintf('Obstacle clearance along the sim (payload radius %.0f mm included):\n', 1000*R_pay);
    for j = 1:numel(S.obs)
        d = arrayfun(@(i) obstacle_sdf(Qg(:, i), S.obs(j), R_pay), 1:size(Qg, 2));
        [dmin, imin] = min(d);
        if dmin < 0, tag = 'COLLISION'; else, tag = 'clear'; end
        fprintf('  obs %d (%s): min %+.1f mm at t = %.2f s  [%s]\n', j, S.obs(j).type, 1000*dmin, t(imin), tag);
    end
end

%% 3. Plots (solid = sim, dots = optimizer nodes)
figure('Name', 'Crane trajectory', 'Color', 'w', 'Position', [60 60 1150 750]);
if isfield(S, 'w_t'), wt_str = sprintf('   (free t_f, w_t = %g)', S.w_t); else, wt_str = ''; end
sgtitle(sprintf('Optimal t_f = %.3f s%s', S.tf, wt_str), 'FontWeight', 'bold');
sub = @(k, rows, scale, lbl, ttl) plot_states(k, t, X(rows, :)*scale, S.t, S.x(rows, :)*scale, lbl, ttl, S.tf);
sub(1, 1:2, 1000,        'mm',  'Trolley x, cable l');
sub(2, 3,   180/pi,      'deg', 'Slew \theta');
sub(3, 4:5, 180/pi,      'deg', 'Swing \phi, \psi');
sub(4, 6:7, 1000,        'mm/s', 'Rates dx, dl');

subplot(2, 3, 5); hold on; grid on;
plot(t, U(1, :)*1000, t, U(2, :)*1000, t, rad2deg(U(3, :)), 'LineWidth', 1.5);
xline(S.tf, 'k-.', sprintf('t_f = %.2f s', S.tf));  xlabel('t (s)');  ylabel('mm/s^2, deg/s^2');  title('Inputs u');
legend('ddx', 'ddl', 'dd\theta', 'Location', 'best');

subplot(2, 3, 6); hold on; grid on; axis equal;
P = [X(1, :).*cos(X(3, :)); X(1, :).*sin(X(3, :))];
plot(1000*P(1, :), 1000*P(2, :), 'Color', [.6 .6 .6], 'LineWidth', 1.2);
plot(1000*Q(1, :), 1000*Q(2, :), 'b', 'LineWidth', 1.5);
plot(1000*Q(1, 1), 1000*Q(2, 1), 'go', 'MarkerFaceColor', 'g');
plot(1000*Qf(1), 1000*Qf(2), 'ks', 'MarkerFaceColor', 'k');
lgd = {'trolley', 'payload', 'start', 'end'};
if has_obs
    for j = 1:numel(S.obs)
        draw_obstacle_top(S.obs(j));
        lgd{end+1} = sprintf('obs %d', j); %#ok<SAGROW>
    end
end
xlabel('X (mm)');  ylabel('Y (mm)');  title('Top view');  legend(lgd, 'Location', 'best');

%% 4. Animation
if animate
    tA = 0 : 1/fps : t(end);
    XA = interp1(t, X.', tA).';
    QA = payload_xyz(XA);  PA = [XA(1, :).*cos(XA(3, :)); XA(1, :).*sin(XA(3, :)); zeros(size(tA))];
    lift = [0; 0; H];
    fa = figure('Name', 'Crane animation', 'Color', 'w', 'WindowStyle', 'normal', ...
                'Position', [1230 60 700 750]);                % own window, not docked in the editor
    ax = axes(fa);  hold(ax, 'on');  grid(ax, 'on');
    axis(ax, 'equal');  view(ax, 40, 25);
    R = 1.2*max(S.x_ub(1), 0.5);
    xlim(ax, [-R R]);  ylim(ax, [-R R]);  zlim(ax, [0 H + 0.1]);  xlabel(ax, 'X (m)');  ylabel(ax, 'Y (m)');  zlabel(ax, 'Z (m)');
    plot3(ax, [0 0], [0 0], [0 H], 'Color', [.3 .3 .3], 'LineWidth', 5);
    if has_obs
        for j = 1:numel(S.obs), draw_obstacle_3d(ax, S.obs(j)); end
    end
    hj = plot3(ax, [0 0], [0 0], [H H], 'Color', [.85 .6 .1], 'LineWidth', 3);
    hc = plot3(ax, [0 0], [0 0], [0 0], 'k-', 'LineWidth', 1.5);
    hp = plot3(ax, 0, 0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
    ht = plot3(ax, 0, 0, 0, 'r-');
    for k = 1:numel(tA)
        if ~isvalid(fa), break; end
        th = XA(3, k);  Pk = PA(:, k) + lift;  Qk = QA(:, k) + lift;
        set(hj, 'XData', [0 S.x_ub(1)*cos(th)], 'YData', [0 S.x_ub(1)*sin(th)], 'ZData', [H H]);
        set(hc, 'XData', [Pk(1) Qk(1)], 'YData', [Pk(2) Qk(2)], 'ZData', [Pk(3) Qk(3)]);
        set(hp, 'XData', Qk(1), 'YData', Qk(2), 'ZData', Qk(3));
        set(ht, 'XData', QA(1, 1:k), 'YData', QA(2, 1:k), 'ZData', QA(3, 1:k) + H);
        if tA(k) > S.tf, phase = 'coast (u = 0)'; else, phase = 'maneuver'; end
        title(ax, sprintf('%s   t = %.2f s', phase, tA(k)));
        drawnow;  pause(1/fps);
    end
end

%% Local functions
function Q = payload_xyz(X)
% Payload position in the ground frame (z relative to the jib, negative down), 3 x N.
x = X(1, :);  l = X(2, :);  th = X(3, :);  phi = X(4, :);  psi = X(5, :);
rx = x + l.*sin(phi).*cos(psi);                 % in the jib frame
ry = l.*sin(psi);
Q  = [rx.*cos(th) - ry.*sin(th);  rx.*sin(th) + ry.*cos(th);  -l.*cos(phi).*cos(psi)];
end

function plot_states(k, t, Y, tn, Yn, lbl, ttl, tf)
subplot(2, 3, k);  hold on;  grid on;
plot(t, Y, 'LineWidth', 1.5);
set(gca, 'ColorOrderIndex', 1);
plot(tn, Yn, '.', 'MarkerSize', 10);
xline(tf, 'k-.', sprintf('t_f = %.2f s', tf));  xlabel('t (s)');  ylabel(lbl);  title(ttl);
end

function draw_obstacle_top(ob)
% Obstacle footprint on the top-view plot (mm).
c = [0.85 0.2 0.2];
switch ob.type
    case 'box'
        lo = ob.center(1:2) - ob.half(1:2);  w = 2*ob.half(1:2);
        rectangle('Position', 1000*[lo.' w.'], 'FaceColor', [c 0.35], 'EdgeColor', c);
        plot(NaN, NaN, 's', 'MarkerFaceColor', c, 'MarkerEdgeColor', c);   % legend proxy
    case 'cylinder'
        a = linspace(0, 2*pi, 60);
        fill(1000*(ob.xy(1) + ob.R*cos(a)), 1000*(ob.xy(2) + ob.R*sin(a)), c, 'FaceAlpha', 0.35, 'EdgeColor', c);
end
end

function draw_obstacle_3d(ax, ob)
% Obstacle as a translucent solid in the animation (m, Z up from the floor).
c = [0.85 0.2 0.2];
switch ob.type
    case 'box'
        lo = ob.center - ob.half;  hi = ob.center + ob.half;
        V = [lo(1) lo(2) lo(3); hi(1) lo(2) lo(3); hi(1) hi(2) lo(3); lo(1) hi(2) lo(3);
             lo(1) lo(2) hi(3); hi(1) lo(2) hi(3); hi(1) hi(2) hi(3); lo(1) hi(2) hi(3)];
        F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
        patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', c, 'FaceAlpha', 0.35, 'EdgeColor', c);
    case 'cylinder'
        [cx, cy, cz] = cylinder(ob.R, 40);
        surf(ax, cx + ob.xy(1), cy + ob.xy(2), ob.z(1) + cz*(ob.z(2) - ob.z(1)), ...
             'FaceColor', c, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        a = linspace(0, 2*pi, 41);
        fill3(ax, ob.xy(1) + ob.R*cos(a), ob.xy(2) + ob.R*sin(a), ob.z(2)*ones(size(a)), c, ...
              'FaceAlpha', 0.35, 'EdgeColor', c);
end
end
