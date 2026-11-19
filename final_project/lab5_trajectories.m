%% lab5_trajectories.m - Input-shaping test trajectories for the Lab 5 tower crane
% Writes Lab 5 continuous-velocity upload files (lab5/*.xlsx: first sheet, headers in row 1:
% time_s, trolley_pct, slew_pct, hoist_pct; one row every 0.05 s, velocity in % of the Table IV
% maximum, ending at zero) that test whether input shaping removes the swing on the real crane.
% Every test starts from the Lab 5 start pose (trolley 550 mm, slew 195 deg, cable 1350 mm;
% "Return to start" puts the crane there) and makes one point-to-point move, unshaped or shaped.
% The lab crane is a double pendulum (hook on the cable, payload 0.381 m below it), so the
% shapers target both of its swing modes. All times sit on the 0.05 s grid. Each file is
% replayed on a linear double-pendulum model (the same ground-frame model as solve_is.m, with a
% hook and a payload) to predict the swing left after the move and the Lab 5 translation
% penalty: lab5/summary.csv, lab5/predicted.png.
% Test plan: run trolley_unshaped first, read the two swing periods off its recorded data,
% put them in T_meas, rerun this script, then run the shaped files.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here);
out_dir = fullfile(here, 'lab5');
if ~isfolder(out_dir), mkdir(out_dir); end

%% Settings (Lab 5 Tables III and IV, and the start pose)
g  = 9.81;
L1 = 1.350;  L2 = 0.381;              % hook cable and payload rigging lengths (m)
m_hook = 0.357;  m_load = 0.050;      % kg
T_meas = [];                          % measured swing periods [T1 T2] (s); [] = use the model's
P0  = [0.550; 1.350; deg2rad(195)];   % start pose [trolley (m); cable (m); slew (rad)]
lo  = [0.500; 0.650; deg2rad(30)];    % position limits, same order
hi  = [0.925; 1.450; deg2rad(350)];
vmax = [72.6e-3; 30e-3; deg2rad(32.4)];   % m/s, m/s, rad/s
amax = [72.6e-3; 30e-3; deg2rad(32.4)];   % m/s^2, m/s^2, rad/s^2
dt   = 0.05;                          % Lab 5 time grid (s)
t_after = 15;                         % residual vibration is recorded this long after the move
% Tests: name, axis (1 trolley, 3 slew), distance (m or rad), top speed (fraction of the
% limit), shaper. The slew tests turn away from the finish (195 -> 235 deg) to keep clear of
% the T obstacle: check that side of the workspace is clear before running them.
% Shaping the slew removes the swing its acceleration starts, but slewing also flings the
% payload outward (centrifugal, growing with speed squared), which no slew shaper cancels;
% the half-speed test shows that part shrinking.
tests = {
    'trolley_unshaped',    1, 0.300,        1,   'none'
    'trolley_ZV_mode1',    1, 0.300,        1,   'ZV1'
    'trolley_ZV_both',     1, 0.300,        1,   'ZV12'
    'slew_unshaped',       3, deg2rad(40),  1,   'none'
    'slew_ZV_both',        3, deg2rad(40),  1,   'ZV12'
    'slew_ZV_both_half',   3, deg2rad(40),  0.5, 'ZV12'
};

%% Swing modes and shapers
% Small swings of a hook (m1) on a cable L1 with a payload (m2) L2 below it, R = m2/m1:
%   w^4 - g(1+R)(1/L1 + 1/L2) w^2 + g^2 (1+R)/(L1 L2) = 0
Rm = m_load/m_hook;  A = (1 + Rm)*(1/L1 + 1/L2);  B = (1 + Rm)/(L1*L2);
w  = sqrt(g/2*(A + [-1, 1]*sqrt(A^2 - 4*B)));
if ~isempty(T_meas), w = 2*pi./T_meas; end
half = round(pi./w/dt)*dt;            % half periods, rounded onto the 0.05 s grid
fprintf('Swing modes: T1 = %.3f s, T2 = %.3f s (ZV impulse spacing %.2f s and %.2f s)\n', 2*pi./w, half);
zv1 = [1/2, 1/2; 0, half(1)];  zv2 = [1/2, 1/2; 0, half(2)];
shapers = struct('none', [1; 0], 'ZV1', zv1, 'ZV12', conv_shapers(zv1, zv2));

%% Build, check, write and predict every test
rows = {};  pred = struct('name', {}, 't', {}, 'd', {}, 'ax', {});
for k = 1:size(tests, 1)
    [name, ax, D, speed, sh_name] = tests{k, :};
    S = crane.move_steps(0, D, speed*vmax(ax), amax(ax));    % velocity steps of the unshaped move
    S(:, 1) = round(S(:, 1)/dt)*dt;                          % on the 0.05 s grid
    R = crane.shape_steps(ax, S, shapers.(sh_name), amax(ax));
    tbl = continuous_table(R, vmax, dt);
    check_table(tbl, P0, lo, hi, vmax, dt);
    writetable(tbl, fullfile(out_dir, [name '.xlsx']), 'Sheet', 1, 'WriteMode', 'replacefile');
    [r, t, d] = predict(tbl, P0, vmax, Rm, L1, L2, g, dt, t_after);
    pred(end+1) = struct('name', name, 't', t, 'd', d, 'ax', ax);       %#ok<SAGROW>
    move = 1000*D;  unit = 'mm';  if ax == 3, move = rad2deg(D); unit = 'deg'; end
    rows(end+1, :) = {name, sh_name, 100*speed, r.t_move, sprintf('%.0f %s', move, unit), ...
                      r.peak, r.integral, 0.02*r.integral};           %#ok<SAGROW>
end
summary = cell2table(rows, 'VariableNames', {'file', 'shaper', 'speed_pct', 'move_time_s', 'move', ...
    'pred_peak_residual_mm', 'pred_residual_mm_s', 'pred_translation_penalty_s'});
disp(summary);
writetable(summary, fullfile(out_dir, 'summary.csv'));

% Predicted payload offset from its final rest point, after each move ends
fig = figure('Name', 'Predicted residual swing', 'Color', 'w', 'Position', [60, 60, 1200, 460]);
for sp = 1:2
    subplot(1, 2, sp); hold on; grid on;
    for p = pred([pred.ax] == 2*sp - 1)
        plot(p.t, 1000*p.d, 'LineWidth', 1.4, 'DisplayName', strrep(p.name, '_', ' '));
    end
    xlabel('Time after the move ends (s)'); ylabel('Payload offset (mm)'); legend('Location', 'northeast');
    if sp == 1, title('Trolley tests (linear double-pendulum model)'); else, title('Slew tests'); end
end
exportgraphics(fig, fullfile(out_dir, 'predicted.png'), 'Resolution', 150);
fprintf('Wrote %d trajectories, summary.csv and predicted.png to %s\n', size(tests, 1), out_dir);

%% Local functions
function sh = conv_shapers(a, b)
% Convolution of two shapers: every pair of impulses, amplitudes multiplied, times added
    A = a(1,:).' * b(1,:);  T = a(2,:).' + b(2,:);
    [t, k] = sort(T(:).');  amp = A(:).';
    sh = [amp(k); t];
end

function tbl = continuous_table(R, vmax, dt)
% Lab 5 continuous file: the average velocity over each 0.05 s step (so the positions match the
% ramp command exactly at every step), as % of the maximum, ending with zero rows
    t_end = max(R(:,2) + R(:,3));
    t = (0:dt:(ceil(t_end/dt) + 2)*dt).';
    P = crane.profile(R, [t; t(end) + dt], zeros(3, 1));
    V = diff(P, 1, 2)/dt;                                     % 3 x numel(t): trolley, hoist, slew
    pct = round(100*V./vmax, 4);
    tbl = table(t, pct(1,:).', pct(3,:).', pct(2,:).', ...
        'VariableNames', {'time_s', 'trolley_pct', 'slew_pct', 'hoist_pct'});
end

function check_table(tbl, P0, lo, hi, vmax, dt)
% The Lab 5 continuous-file rules, and the position limits along the way
    t = tbl.time_s;
    V = [tbl.trolley_pct, tbl.hoist_pct, tbl.slew_pct].'/100 .* vmax;
    assert(t(1) == 0 && all(abs(diff(t) - dt) < 1e-9) && t(end) <= 80, 'rows must run 0, 0.05, ... up to 80 s');
    assert(all(abs([tbl.trolley_pct; tbl.slew_pct; tbl.hoist_pct]) <= 100), 'velocities must be within 100 %');
    assert(all(V(:, end) == 0), 'must end with zero velocity');
    P = P0 + cumsum(V*dt, 2);
    assert(all(P >= lo - 1e-9 & P <= hi + 1e-9, 'all'), 'position limit exceeded');
end

function [r, t_out, d_out] = predict(tbl, P0, vmax, Rm, L1, L2, g, dt, t_after)
% Replay an upload table on a linear double-pendulum model in the ground frame: hook z1 on L1
% below the trolley p = x e^{i theta}, payload z2 on L2 below the hook (complex numbers),
%   z1'' = (1+R) g/L1 (p - z1) + R g/L2 (z2 - z1),   z2'' = g/L2 (z1 - z2)
% Each row's velocity is held for 0.05 s. Returns the move time and the payload's offset from
% its final rest point afterwards: peak (mm) and integral (mm*s) over t_after seconds.
    t = tbl.time_s;
    V = [tbl.trolley_pct, tbl.hoist_pct, tbl.slew_pct].'/100 .* vmax;
    t_move = t(find(any(V ~= 0, 1), 1, 'last')) + dt;
    edges = [t; t(end) + dt];
    Pe = P0 + [zeros(3, 1), cumsum(V*dt, 2)];                 % exact positions at the step edges
    pos = @(tt) interp1(edges, Pe.', min(tt(:), edges(end))).';
    T_end = t_move + t_after;  h = 1e-3;  N = round(T_end/h);  h = T_end/N;
    tt = linspace(0, T_end, 2*N + 1);                         % RK4 grid with half steps
    P = pos(tt);  p = P(1,:) .* exp(1i*P(3,:));
    c1 = (1 + Rm)*g/L1;  c2 = Rm*g/L2;  c3 = g/L2;
    f = @(pp, s) [s(2); c1*(pp - s(1)) + c2*(s(3) - s(1)); s(4); c3*(s(1) - s(3))];
    s = [p(1); 0; p(1); 0];  z2 = complex(zeros(1, N + 1));  z2(1) = s(3);
    for k = 1:N
        p0 = p(2*k-1);  pm = p(2*k);  p1 = p(2*k+1);
        k1 = f(p0, s);  k2 = f(pm, s + h/2*k1);  k3 = f(pm, s + h/2*k2);  k4 = f(p1, s + h*k3);
        s = s + h/6*(k1 + 2*k2 + 2*k3 + k4);
        z2(k+1) = s(3);
    end
    ts = tt(1:2:end);  after = ts >= t_move - 1e-9;
    d = abs(z2(after) - p(end));                              % offset from the final rest point
    t_out = ts(after) - t_move;  d_out = d;
    r.t_move = t_move;
    r.peak = 1000*max(d);
    r.integral = 1000*trapz(t_out, d);
end
