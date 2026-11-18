function M = crane_metrics(T, cfg)
% Limit usage, pin clearance and optimality checks for a trajectory from crane_trajectory.
% Ratios are against the PLANNING limits (hardware x limit_scale); M.admissible is true
% only if nothing is exceeded and the payload never touches the clearance circle.
X = T.X;  G = T.G;
M.t_f         = T.t(end);
M.H_spread    = max(T.H) - min(T.H);                  % PMP: H is constant along an optimum (autonomous)
M.x_range     = [min(X(1,:)), max(X(1,:))];
M.l_range     = [min(X(2,:)), max(X(2,:))];
M.th_range_deg = rad2deg([min(X(3,:)), max(X(3,:))]) + cfg.pin_slew_deg;   % real slew angles
M.v_ratio     = [max(abs(X(6,:)))/cfg.dxmax, max(abs(X(7,:)))/cfg.dlmax, max(abs(X(8,:)))/cfg.dthmax];
M.u_ratio     = (max(abs(T.U), [], 2) ./ [cfg.ddxmax; cfg.ddlmax; cfg.ddthmax]).';
M.d_pin_min   = min(G.d_pin);
M.d_cable_min = min(G.d_cable);
M.reach_over  = max(G.r_rad) - cfg.xmax;              % > 0: the payload went where the trolley can't
M.sweep_deg   = rad2deg(X(11, end));
M.swing_max_deg = rad2deg(max(abs(X(4:5, :)), [], 2)).';

tol = cfg.tol;
M.why = {};
if M.x_range(1) < cfg.xmin*(1-tol) || M.x_range(2) > cfg.xmax*(1+tol), M.why{end+1} = 'x'; end
if M.l_range(1) < cfg.lmin*(1-tol) || M.l_range(2) > cfg.lmax*(1+tol), M.why{end+1} = 'x'; end
if any(M.v_ratio > 1 + tol),                                          M.why{end+1} = 'velocity'; end
if M.d_pin_min < (cfg.r_p + cfg.c_clr)*(1 - tol),                     M.why{end+1} = 'payload-pin'; end
if M.d_cable_min < cfg.r_p,                                           M.why{end+1} = 'cable-pin'; end
if abs(X(11, end) - cfg.sweep) > 1e-2,                                M.why{end+1} = 'wrong side of pin'; end
M.why = unique(M.why, 'stable');
M.admissible = isempty(M.why);
end
