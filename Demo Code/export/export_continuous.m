function cmd = export_continuous(S, file, tune)
% Upload file for the crane GUI: time_s, trolley_pct, slew_pct, hoist_pct on a 0.05 s grid (% of max velocity).
% Row k holds the velocity from t_k to t_k+1. It is the AVERAGE planned velocity over that
% interval, (q(t_k+1) - q(t_k))/dt, so integrating the commands lands exactly on the planned
% path at every grid point (no zero-order-hold drift). The last row is all zeros.
%
% tune.axis_gain = [trolley slew hoist]: measured (actual move / commanded move). The
% commands are divided by it, so a drive that overshoots by 2 % (gain 1.02) is sent 2 % less.
dt = 0.05;
t  = (0:dt:ceil(S.t_f/dt - 1e-9)*dt).';
Z  = crane_eval(S, t.');                                   % holds the final rest pose past t_f
q  = Z(1:3, :).';                                          % [x l θ] in m, m, rad

v   = [diff(q) ./ dt; 0 0 0];                              % average velocity per row, last row = stop
hw  = S.cfg.hw;
vmax = [hw.dxmax, hw.dthmax, hw.dlmax];                    % trolley, slew, hoist (SI)
pct = 100 * v(:, [1 3 2]) ./ vmax ./ tune.axis_gain;       % columns: trolley, slew (+ = slew angle up), hoist (+ = lower)

over = max(abs(pct), [], 1);
if any(over > 100)
    warning('Commands exceed 100 %% (trolley %.1f, slew %.1f, hoist %.1f): clipped. Lower limit_scale.', over);
    pct = max(min(pct, 100), -100);
end

T = array2table([t, round(pct, 4)], 'VariableNames', {'time_s', 'trolley_pct', 'slew_pct', 'hoist_pct'});
if isfile(file), delete(file); end                         % writetable would keep stale rows of an old file
writetable(T, file, 'Sheet', 1);

cmd.file = file;  cmd.t = t;  cmd.pct = pct;
cmd.v = pct/100 .* vmax .* tune.axis_gain;                 % what the crane will actually do (SI, [x θ l])
cmd.peak_pct = over;
end
