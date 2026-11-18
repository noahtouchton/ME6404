function plot_preview(S, R)
% Plan vs replay of the upload file: top view, drive velocities, pin gap, residual swing.
c   = S.cfg;  hw = c.hw;
T   = crane_trajectory(S, linspace(0, S.t_f, 800));
rot = deg2rad(c.pin_slew_deg);
Rz  = [cos(rot) -sin(rot); sin(rot) cos(rot)];            % model frame -> real slew angles
w   = @(XY) 1000 * Rz * XY;                               % to mm, real orientation

f = findobj('Type', 'figure', 'Tag', 'lab_preview');
if isempty(f), f = figure('Tag', 'lab_preview', 'Name', 'Lab preview', 'Color', 'w'); end
clf(f);  tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
title(tl, sprintf('t_f = %.2f s   |   replay of the upload file: residual swing %.2f mm', S.t_f, R.residual_peak_mm));

% --- Top view ---
ax = nexttile(tl);  hold(ax, 'on');  axis(ax, 'equal');  grid(ax, 'on');
ang = linspace(0, 2*pi, 100);
th  = T.X(3, :) + rot;  arc = linspace(min(th) - 0.15, max(th) + 0.15, 100);
plot(ax, 1000*hw.xmax*cos(arc), 1000*hw.xmax*sin(arc), ':', 'Color', [.6 .6 .6]);
pin = w([c.pin_x; c.pin_y]);
fill(ax, pin(1) + 1000*c.r_p*cos(ang), pin(2) + 1000*c.r_p*sin(ang), [.3 .3 .3], 'EdgeColor', 'none');
plot(ax, pin(1) + 1000*(c.r_p + c.c_clr)*cos(ang), pin(2) + 1000*(c.r_p + c.c_clr)*sin(ang), 'r--');
Pp = w(T.G.Qp(1:2, :));  Pt = w(T.G.P(1:2, :));  Pr = w(R.G.Qp(1:2, :));
plot(ax, Pt(1,:), Pt(2,:), 'Color', [.6 .6 .6], 'LineWidth', 1);
plot(ax, Pp(1,:), Pp(2,:), 'b', 'LineWidth', 1.5);
plot(ax, Pr(1,:), Pr(2,:), 'm--', 'LineWidth', 1);
plot(ax, Pp(1,1), Pp(2,1), 'go', 'MarkerFaceColor', 'g');  plot(ax, Pp(1,end), Pp(2,end), 'ks', 'MarkerFaceColor', 'k');
legend(ax, {'trolley limit', 'pin', 'clearance', 'trolley', 'payload (plan)', 'payload (replay)', 'start', 'finish'}, 'Location', 'best');
xlabel(ax, 'mm');  ylabel(ax, 'mm');  title(ax, 'Top view');

% --- Drive velocities, % of the hardware limit ---
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
vm = [hw.dxmax; hw.dthmax; hw.dlmax];  col = lines(3);
plot(ax, T.t, 100*T.X([6 8 7], :) ./ vm, 'LineWidth', 1.5);
set(ax, 'ColorOrderIndex', 1);
plot(ax, R.t, 100*R.X([6 8 7], :) ./ vm, '--');
legend(ax, {'trolley', 'slew', 'hoist'}, 'Location', 'best');  xlabel(ax, 't (s)');  ylabel(ax, '% of max velocity');
title(ax, 'Velocities (solid plan, dashed replay)');  ax.ColorOrder = col;

% --- Pin gap ---
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
plot(ax, T.t, 1000*(T.G.d_pin - c.r_p), 'b', R.t, 1000*(R.G.d_pin - c.r_p), 'm--');
plot(ax, T.t, 1000*(T.G.d_cable - c.r_p), 'Color', [.4 .7 .4]);
yline(ax, 1000*c.c_clr, 'r--', 'clearance');  yline(ax, 0, 'k');
ylim(ax, [-5, max(100, 3000*c.c_clr)]);  xlabel(ax, 't (s)');  ylabel(ax, 'gap to pin surface (mm)');
legend(ax, {'payload (plan)', 'payload (replay)', 'cable (plan)'}, 'Location', 'best');  title(ax, 'Pin clearance');

% --- Residual swing after the stop ---
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
Xt = c.X_end;  Xt(11) = 0;  Qt = crane_geometry(Xt, c).Qp(1:2);
plot(ax, R.t, R.offset_mm, 'm');
xline(ax, R.t_cmd, 'k--', 'motion ends');
xlabel(ax, 't (s)');  ylabel(ax, 'payload distance from target (mm)');
title(ax, sprintf('Replay: residual swing %.2f mm after the stop', R.residual_peak_mm));
end
