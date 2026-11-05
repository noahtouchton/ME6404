function M = crane_metrics(T, cfg)
% Limit usage, pin clearance and optimality checks for a trajectory from crane_trajectory.
% M.admissible is true only if the maneuver respects every hardware limit and never
% touches the pin; M.why lists what failed.
X = T.X;  G = T.G;
M.t_f         = T.t(end);
M.H_spread    = max(T.H) - min(T.H);
M.x_ratio     = max(X(1,:)) / cfg.xmax;
M.v_ratio     = [max(abs(X(6,:)))/cfg.dxmax, max(abs(X(7,:)))/cfg.dlmax, max(abs(X(8,:)))/cfg.dthmax];
M.u_ratio     = (max(abs(T.U), [], 2) ./ [cfg.ddxmax; cfg.ddlmax; cfg.ddthmax]).';
M.d_pin_min   = min(G.d_pin);
M.d_cable_min = min(G.d_cable);
M.reach_over  = max(G.r_rad) - cfg.xmax;
M.windings    = X(11, end) / (2*pi);
M.end_swing   = rad2deg(X(4:5, end)).';

tol = 1e-3;                                         % 0.1 % slack for the soft walls
M.why = {};
if M.x_ratio > 1 + tol,                          M.why{end+1} = 'x'; end
if any(M.v_ratio > 1 + tol),                     M.why{end+1} = 'velocity'; end
if M.d_pin_min < (cfg.r_p + cfg.c_clr)*(1 - tol), M.why{end+1} = 'payload-pin'; end
if M.d_cable_min < cfg.r_p,                      M.why{end+1} = 'cable-pin'; end
if abs(M.windings - cfg.kdir) > 1e-3,            M.why{end+1} = 'windings'; end
M.admissible = isempty(M.why);
end
