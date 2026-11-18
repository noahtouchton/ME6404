function [guess, t_f] = initial_guess(cfg)
% Cold-start guess: the payload sweeps cfg.sweep radians around the pin center (pulled in to
% the clearance circle halfway), hanging straight down, at rest at both ends. Only the
% topology has to be right: Newton fixes the rest.
%   * The trolley radius is clipped to its range, so the guess never asks for the impossible.
%   * Timing comes from the velocity limits segment by segment (each piece of path takes as
%     long as its slowest axis needs at guess_slack x below the limit), with smooth ramps at
%     the ends, so the horizon is as short as the path allows. Long horizons make the
%     costates stiff and the solve fail.
%   * Costates get a tiny random value: exactly zero makes λᵀ∂f/∂x vanish (singular Newton matrix).
N = 400;
s = linspace(0, 1, N);
s = 0.5 - 0.5*cos(pi*s);                                   % denser samples near the ends

% Polar path around the pin center
pin  = [cfg.pin_x; cfg.pin_y];
Q0   = cfg.X_start(1) * [cos(cfg.X_start(3)); sin(cfg.X_start(3))];
Qf   = cfg.X_end(1)   * [cos(cfg.X_end(3));   sin(cfg.X_end(3))];
a0   = atan2(Q0(2) - pin(2), Q0(1) - pin(1));
rho0 = norm(Q0 - pin);   rhof = norm(Qf - pin);
rhoc = cfg.r_p + cfg.c_clr + cfg.guess_margin;             % closest approach
rho_lin = rho0 + (rhof - rho0)*s;
rho  = rho_lin - max(rho_lin - rhoc, 0) .* sin(pi*s).^2;   % pulled in to rhoc at mid-path
Q    = pin + rho .* [cos(a0 + cfg.sweep*s); sin(a0 + cfg.sweep*s)];

x   = min(max(vecnorm(Q), cfg.xmin + 0.005), cfg.xmax - 0.005);
th  = unwrap(atan2(Q(2,:), Q(1,:)));
l   = cfg.X_start(2) + (cfg.X_end(2) - cfg.X_start(2))*s;
q   = [x; l; th];

% Time along the path. p = "seconds at full speed" (each segment priced by its slowest
% axis), then a trapezoidal profile in p: dp/dt ≤ 1/slack, d²p/dt² ≤ a_n/slack, where
% a_n = min(a_max/v_max) over the axes (= 1/s for this crane).
vmax = [cfg.dxmax; cfg.dlmax; cfg.dthmax];
amax = [cfg.ddxmax; cfg.ddlmax; cfg.ddthmax];
p    = [0, cumsum(max(abs(diff(q, 1, 2)) ./ vmax, [], 1))];
P    = p(end);
v    = 1/cfg.guess_slack;  a = min(amax ./ vmax)/cfg.guess_slack;
pr   = min(v^2/(2*a), P/2);                                % path used by each ramp
vp   = sqrt(2*a*pr);                                       % peak speed (< v if no cruise)
t_f  = 2*vp/a + (P - 2*pr)/v;
t    = zeros(1, N);
up   = p <= pr;  dn = p >= P - pr;  cr = ~up & ~dn;
t(up) = sqrt(2*p(up)/a);
t(cr) = vp/a + (p(cr) - pr)/v;
t(dn) = t_f - sqrt(2*max(P - p(dn), 0)/a);
[t, keep] = unique(t);                                     % drop repeated samples (no motion)
q = q(:, keep);  s = s(:, keep);  N = numel(t);
x = q(1, :);  l = q(2, :);  th = q(3, :);

Y = zeros(22, N);
Y(1:3, :) = q;
Y(6:8, :) = [gradient(x, t); gradient(l, t); gradient(th, t)];
Y(6:8, [1 end]) = 0;
Y(11, :)  = cfg.sweep*s;
rng(cfg.seed);
Y(12:22, :) = 1e-3 * randn(11, N);

guess = bvpinit(t/t_f, zeros(22, 1));
guess.y = Y;
end
