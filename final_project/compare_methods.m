%% compare_methods.m - Optimal control vs input shaping on the same loop around the pin
% Replays each saved solution's actuator command open loop on the full nonlinear model
% (optimal control: U*(t) from the bvp4c costates, so solver tolerance shows up too), lets the
% payload coast 5 s with the actuators stopped, and repeats the replay with the wrong cable
% length. Prints a table and writes results/comparison.csv, comparison.mat and figures/*.png.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'model', 'generated'));
t_start = tic;

%% Settings
files    = {'oc_fixed', 'oc_min_time', 'is_classic', 'is_cancel', 'is_robust'};   % results/<name>.mat
l_true   = 0.80:0.02:1.20;    % real cable length (m); every design assumed l0 = 1.00 m
t_coast  = 5;                 % s with the actuators stopped after the maneuver
ok_swing = 0.005;             % leftover swing (m) that still counts as at rest
fig_dir  = fullfile(fileparts(crane.file('comparison.csv')), 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
stopped  = struct('U', @(t) zeros(3, numel(t)), 'breaks', [], 'pwc', false);

%% 1. Replay each solution: as designed, coasting, and with the wrong cable length
res = struct('S', {}, 'T', {}, 'Tc', {}, 'M', {}, 'coast_swing', {}, 'sweep', {}, 'band', {});
for f = files
    p = crane.file([f{1} '.mat']);
    if ~isfile(p), fprintf('skipping %s: run solve_oc.m / solve_is.m first\n', f{1}); continue; end
    S = load(p);  c = crane.command(S);  cfg = S.cfg;
    t  = linspace(0, S.t_f, 3000);
    T  = crane.simulate(c, t, cfg);
    Tc = crane.simulate(stopped, S.t_f + linspace(0, t_coast, 400), cfg, T.X(:, end));
    coast_swing = max(vecnorm(Tc.G.Qp(1:2,:) - Tc.G.P(1:2,:)));   % payload offset from the stopped trolley
    sw = struct('resid', nan(size(l_true)), 'pin', nan(size(l_true)), 'adm', false(size(l_true)));
    for k = 1:numel(l_true)
        X0 = [cfg.X0; 0];  X0(2) = l_true(k);
        Ml = crane.metrics(crane.simulate(c, t, cfg, X0), cfg);
        sw.resid(k) = Ml.residual;
        sw.pin(k)   = Ml.d_pin_min - cfg.r_p - cfg.c_clr;
        sw.adm(k)   = Ml.admissible;
    end
    good = sw.adm & sw.resid <= ok_swing;                 % robust band: contiguous around l0
    [~, k0] = min(abs(l_true - cfg.l0));  lo = k0;  hi = k0;  band = [NaN, NaN];
    if good(k0)
        while lo > 1 && good(lo-1), lo = lo - 1; end
        while hi < numel(l_true) && good(hi+1), hi = hi + 1; end
        band = l_true([lo, hi]);
    end
    res(end+1) = struct('S', S, 'T', T, 'Tc', Tc, 'M', crane.metrics(T, cfg), ...
                        'coast_swing', coast_swing, 'sweep', sw, 'band', band);  %#ok<SAGROW>
    fprintf('%-26s replayed (%.0f s so far)\n', S.label, toc(t_start));
end
n = numel(res);
if n == 0, error('No solutions in results/ yet: run solve_oc.m and solve_is.m first.'); end

%% 2. Table
lab = arrayfun(@(r) r.S.label, res, 'UniformOutput', false);
fprintf('\n%-26s %6s %6s %6s %6s %16s %16s %7s %11s %8s\n', 'method', 't_f', 'reach', 'pin', 'cable', ...
    '|v|/vmax x l th', '|a|/amax x l th', 'left', 'band (l, m)', 'compute');
fprintf('%-26s %6s %6s %6s %6s %16s %16s %7s %11s %8s\n', '', '(s)', '(mm)', '(mm)', '(mm)', '', '', '(mm)', '', '(s)');
rows = cell(n, 1);
for i = 1:n
    S = res(i).S;  M = res(i).M;  cfg = S.cfg;
    pin = 1000*(M.d_pin_min - cfg.r_p - cfg.c_clr);  cable = 1000*(M.d_cable_min - cfg.r_p);
    adm = '';  if ~M.admissible, adm = ['  NOT ADMISSIBLE: ' strjoin(M.why, ', ')]; end
    fprintf('%-26s %6.2f %6.1f %6.1f %6.1f %16s %16s %7.2f %11s %8.1f%s\n', lab{i}, S.t_f, 1000*M.reach_over, ...
        pin, cable, sprintf('%.2f %.2f %.2f', M.v_ratio), sprintf('%.2f %.2f %.2f', M.u_ratio), ...
        1000*res(i).coast_swing, sprintf('%.2f-%.2f', res(i).band), S.compute_time, adm);
    rows{i} = {lab{i}, S.t_f, 1000*M.reach_over, pin, cable, M.v_ratio(1), M.v_ratio(2), M.v_ratio(3), ...
               M.u_ratio(1), M.u_ratio(2), M.u_ratio(3), 1000*res(i).coast_swing, res(i).band(1), res(i).band(2), ...
               S.compute_time, M.admissible};
end
fprintf(['reach: payload past xmax | pin, cable: margin beyond the required clearance | left: largest payload\n' ...
         'offset from the stopped trolley in the %g s coast | band: cable lengths where the replay is admissible\n' ...
         'and leaves < %.0f mm of swing (every design assumed l0 = %.2f m)\n'], t_coast, 1000*ok_swing, res(1).S.cfg.l0);
tbl = cell2table(vertcat(rows{:}), 'VariableNames', {'method', 't_f_s', 'reach_past_xmax_mm', 'pin_margin_mm', ...
    'cable_margin_mm', 'v_ratio_x', 'v_ratio_l', 'v_ratio_th', 'a_ratio_x', 'a_ratio_l', 'a_ratio_th', ...
    'swing_left_mm', 'band_l_min_m', 'band_l_max_m', 'compute_s', 'admissible'});
writetable(tbl, crane.file('comparison.csv'));
save(crane.file('comparison.mat'), 'tbl', 'l_true', 'res', '-v7.3');

%% 3. Figures (optimal control solid, input shaping dashed)
cols = lines(n);
sty  = repmat({'-'}, 1, n);
sty(arrayfun(@(r) strcmp(r.S.kind, 'commands'), res)) = {'--'};
cfg  = res(1).S.cfg;
ang  = linspace(0, 2*pi, 200);

% Payload paths, top-down: the whole loop and the pass around the pin
fig = figure('Name', 'Payload paths', 'Color', 'w', 'Position', [50, 50, 1200, 560]);
for sp = 1:2
    subplot(1, 2, sp); hold on; grid on; axis equal;
    fill((cfg.R_p + cfg.r_p*cos(ang))*1000, cfg.r_p*sin(ang)*1000, [0.2 0.6 0.2], 'EdgeColor', 'none', 'DisplayName', 'Pin');
    plot((cfg.R_p + (cfg.r_p + cfg.c_clr)*cos(ang))*1000, (cfg.r_p + cfg.c_clr)*sin(ang)*1000, 'g--', 'DisplayName', 'Clearance');
    plot(cfg.xmax*1000*cos(ang), cfg.xmax*1000*sin(ang), 'k--', 'DisplayName', 'x_{max}');
    for i = 1:n
        Q = res(i).T.G.Qp;
        plot(Q(1,:)*1000, Q(2,:)*1000, sty{i}, 'Color', cols(i,:), 'LineWidth', 1.4, 'DisplayName', lab{i});
    end
    xlabel('X (mm)'); ylabel('Y (mm)');
    if sp == 1, title('Payload paths (top-down)');
    else, xlim([820 1000]); ylim([-110 110]); title('Passing the pin'); legend('Location', 'eastoutside'); end
end
exportgraphics(fig, fullfile(fig_dir, 'paths.png'), 'Resolution', 150);

% Swing: payload offset from the trolley, maneuver + coast
fig = figure('Name', 'Swing offset', 'Color', 'w', 'Position', [80, 80, 1000, 450]); hold on; grid on;
for i = 1:n
    t   = [res(i).T.t, res(i).Tc.t];
    off = [vecnorm(res(i).T.G.Qp(1:2,:) - res(i).T.G.P(1:2,:)), vecnorm(res(i).Tc.G.Qp(1:2,:) - res(i).Tc.G.P(1:2,:))];
    plot(t, off*1000, sty{i}, 'Color', cols(i,:), 'LineWidth', 1.3, 'DisplayName', lab{i});
    xline(res(i).S.t_f, ':', 'Color', cols(i,:), 'HandleVisibility', 'off');
end
xlabel('Time (s)'); ylabel('Payload offset from trolley (mm)');
title('Swing during the loop and after it (dotted: end of each maneuver)'); legend('Location', 'northeast');
exportgraphics(fig, fullfile(fig_dir, 'swing.png'), 'Resolution', 150);

% Robustness to cable-length error
fig = figure('Name', 'Robustness', 'Color', 'w', 'Position', [110, 110, 1200, 460]);
subplot(1, 2, 1); hold on; grid on;
for i = 1:n
    plot(sqrt(cfg.l0 ./ l_true), res(i).sweep.resid*1000, sty{i}, 'Marker', '.', 'Color', cols(i,:), 'LineWidth', 1.3, 'DisplayName', lab{i});
end
yline(ok_swing*1000, 'k:', 'HandleVisibility', 'off');
xlabel('\omega_{actual} / \omega_{design}'); ylabel('Swing left at t_f (mm)');
title('Leftover swing vs frequency error'); legend('Location', 'north');
subplot(1, 2, 2); hold on; grid on;
for i = 1:n
    plot(l_true, res(i).sweep.pin*1000, sty{i}, 'Marker', '.', 'Color', cols(i,:), 'LineWidth', 1.3, 'DisplayName', lab{i});
end
yline(0, 'k-', 'HandleVisibility', 'off');
xlabel('Actual cable length (m)'); ylabel('Payload-pin margin (mm)');
title('Pin margin vs cable length (below 0 = clips the clearance)');
exportgraphics(fig, fullfile(fig_dir, 'robustness.png'), 'Resolution', 150);

% Actuator velocities over time
fig = figure('Name', 'Actuator usage', 'Color', 'w', 'Position', [140, 140, 1000, 700]);
names = {'trolley  dx / dx_{max}', 'hoist  dl / dl_{max}', 'slew  d\theta / d\theta_{max}'};
vmax  = [cfg.dxmax; cfg.dlmax; cfg.dthmax];
for a = 1:3
    subplot(3, 1, a); hold on; grid on;
    for i = 1:n
        plot(res(i).T.t, res(i).T.X(5 + a, :) / vmax(a), sty{i}, 'Color', cols(i,:), 'LineWidth', 1.2, 'DisplayName', lab{i});
    end
    ylabel(names{a}); ylim([-1.1 1.1]);
    if a == 1, title('Actuator velocities (fraction of the limit)'); legend('Location', 'eastoutside'); end
end
xlabel('Time (s)');
exportgraphics(fig, fullfile(fig_dir, 'actuators.png'), 'Resolution', 150);
fprintf('\nSaved %s, comparison.mat and figures in %s (%.0f s)\n', crane.file('comparison.csv'), fig_dir, toc(t_start));
