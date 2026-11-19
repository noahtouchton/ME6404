%% trajectories.m - The input-shaped loops as files the tower crane can run
% Writes each saved input-shaped loop (results/is_*.mat, made by solve_is.m) in the crane's
% continuous-velocity upload format: trajectories/<name>.xlsx, first sheet, headers time_s,
% trolley_pct, slew_pct, hoist_pct, one row every 0.05 s, velocity in % of the maximum, ending at
% zero. Each row is the average velocity over its step, so the crane ends up exactly where the
% design does. Each file is then replayed on the full model as a check: trajectories/summary.csv
% and the projected payload paths, trajectories/projected_paths.png.
% The velocities are relative: start the crane at rest at the task's start (trolley 0.70 m,
% cable 1.00 m) with the jib pointing at the pin (0.90 m out along that line). The loop slews to
% one side first, so leave room for its sweep both ways.
% The designs assume a 1.00 m pendulum. If the real payload swings with another period T, set
% l0 = g*(T/(2*pi))^2 in crane.task, rerun solve_is, then this script.
% To run the optimal-control loops on the crane as well, add 'oc_fixed' and 'oc_min_time'.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'model', 'generated'));
out_dir = fullfile(here, 'trajectories');
if ~isfolder(out_dir), mkdir(out_dir); end
loops = {'is_classic', 'is_cancel', 'is_robust'};      % results/<name>.mat to export
dt = 0.05;                                              % upload time step (s)

rows = {};  paths = {};
for f = loops
    S = load(crane.file([f{1} '.mat']));  cfg = S.cfg;
    vmax = [cfg.dxmax; cfg.dlmax; cfg.dthmax];

    % The file: average velocity over each 0.05 s step, as % of the maximum
    t = (0:dt:(ceil(S.t_f/dt) + 2)*dt).';               % rows, ending with zeros after the loop
    P = positions(S, [t; t(end) + dt]);                 % [trolley; cable; slew] at the step edges
    pct = round(100*diff(P, 1, 2)/dt ./ vmax, 4);
    assert(all(abs(pct) <= 100, 'all') && all(pct(:, end) == 0), '%s: speed over 100 %% or not ending at zero', f{1});
    tbl = table(t, pct(1,:).', pct(3,:).', pct(2,:).', ...
        'VariableNames', {'time_s', 'trolley_pct', 'slew_pct', 'hoist_pct'});
    writetable(tbl, fullfile(out_dir, [f{1} '.xlsx']), 'Sheet', 1, 'WriteMode', 'replacefile');

    % Check: replay the file on the full model, the velocity moving linearly from row to row
    % (each row's value at the middle of its step), then 5 s with the motors stopped
    tk = [0; t + dt/2];  vk = [zeros(3, 1), pct/100 .* vmax];
    R = zeros(0, 4);
    for k = 1:numel(tk) - 1
        a = (vk(:, k+1) - vk(:, k)) / (tk(k+1) - tk(k));
        for ax = find(a ~= 0).'
            R(end+1, :) = [ax, tk(k), tk(k+1) - tk(k), a(ax)];   %#ok<AGROW>
        end
    end
    T = crane.simulate(crane.command(struct('ramps', R)), linspace(0, tk(end) + 5, 3000), cfg);
    M = crane.metrics(T, cfg);
    t_move = t(find(any(pct ~= 0, 1), 1, 'last')) + dt;
    rows(end+1, :) = {f{1}, S.label, t_move, 1000*M.reach_over, ...
        1000*(M.d_pin_min - cfg.r_p - cfg.c_clr), 1000*M.residual, M.admissible};   %#ok<SAGROW>
    paths{end+1} = struct('label', S.label, 'T', T, 't_move', t_move, 'left', M.residual);   %#ok<SAGROW>
end
summary = cell2table(rows, 'VariableNames', {'file', 'loop', 'duration_s', 'reach_past_xmax_mm', ...
    'pin_margin_mm', 'swing_left_mm', 'admissible'});
disp(summary);
writetable(summary, fullfile(out_dir, 'summary.csv'));

% Projected payload paths seen from above: solid while the file runs, dotted for the 5 s after
fig = figure('Name', 'Projected paths', 'Color', 'w', 'Position', [60, 60, 1300, 640]);
cols = lines(numel(paths));  ang = linspace(0, 2*pi, 200);  cfg = crane.task();
for sp = 1:2
    subplot(1, 2, sp); hold on; grid on; axis equal;
    fill((cfg.R_p + cfg.r_p*cos(ang))*1000, cfg.r_p*sin(ang)*1000, [0.2 0.6 0.2], 'EdgeColor', 'none', 'DisplayName', 'pin');
    plot((cfg.R_p + (cfg.r_p + cfg.c_clr)*cos(ang))*1000, (cfg.r_p + cfg.c_clr)*sin(ang)*1000, 'g--', 'DisplayName', 'pin clearance');
    plot(cfg.xmax*1000*cos(ang), cfg.xmax*1000*sin(ang), 'k--', 'DisplayName', 'trolley limit');
    for k = 1:numel(paths)
        p = paths{k};  Q = 1000*p.T.G.Qp;  mv = p.T.t <= p.t_move;
        plot(Q(1,mv), Q(2,mv), '-', 'Color', cols(k,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s: %.2f s, %.1f mm left', p.label, p.t_move, 1000*p.left));
        plot(Q(1,~mv), Q(2,~mv), ':', 'Color', cols(k,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    plot(1000*cfg.x0, 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'start and end');
    xlabel('X (mm)'); ylabel('Y (mm)');
    if sp == 1
        title('Projected payload paths (full model, replaying each file)');
        legend('Location', 'southoutside');
    else
        xlim([820 1000]); ylim([-110 110]); title('Passing the pin');
    end
end
exportgraphics(fig, fullfile(out_dir, 'projected_paths.png'), 'Resolution', 150);
fprintf('Wrote %d trajectories, summary.csv and projected_paths.png to %s\n', numel(loops), out_dir);

function P = positions(S, te)
% Trolley, cable and slew positions of a saved loop at times te, held after the loop ends
    if strcmp(S.kind, 'commands')
        P = crane.profile(S.cmd.ramps, te, [S.cfg.x0; S.cfg.l0; 0]);
    else
        T = crane.trajectory(S, min(te(:).', S.t_f));
        P = T.X(1:3, :);
    end
end
