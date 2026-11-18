function S = target_sweep(cfg)
% Total angle (rad) the payload sweeps around the pin center, start -> finish.
%   cfg.passes    times the payload goes behind the pin (crosses the far side: the ray
%                 pointing away from the tower). 0 = stay on the tower side.
%   cfg.direction 'auto' (shortest such path), '+' or '-' (sense of the sweep around the
%                 pin, + = the way the slew angle increases).
a0 = payload_angle(cfg, cfg.X_start);
af = payload_angle(cfg, cfg.X_end);
best = inf;
for k = -4:4
    Sk = af - a0 + 2*pi*k;
    lo = min(a0, a0 + Sk);  hi = max(a0, a0 + Sk);
    n  = sum(2*pi*(-4:4) > lo & 2*pi*(-4:4) < hi);   % crossings of angle 0 (mod 2π)
    dir_ok = strcmp(cfg.direction, 'auto') || sign(Sk) == (1 - 2*strcmp(cfg.direction, '-'));
    if n == cfg.passes && dir_ok && abs(Sk) < abs(best), best = Sk; end
end
assert(isfinite(best), 'No path with %d passes behind the pin in direction ''%s''', cfg.passes, cfg.direction);
S = best;
end

function a = payload_angle(cfg, X)
% Angle of the hanging payload around the pin center, in (-π, π]; π = tower side.
Q = X(1) * [cos(X(3)); sin(X(3))];
a = atan2(Q(2) - cfg.pin_y, Q(1) - cfg.pin_x);
end
