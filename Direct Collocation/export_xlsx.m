function export_xlsx(file, t, u, u_max)
% Write the optimal inputs to an .xlsx with the header
%   Continuous | time_s | trolley_pct | slew_pct | hoist_pct
% One row per time step of the solution (all nodes, no resampling). The first column is
% left blank under "Continuous". The accelerations are given as % of their max:
%   trolley_pct = 100*ddx/ddx_max,  slew_pct = 100*ddth/ddth_max,  hoist_pct = 100*ddl/ddl_max
%
% t     : 1 x n (or n x 1) times (s)
% u     : 3 x n inputs [ddx; ddl; ddth] (SI units)
% u_max : 3 x 1 max accelerations, same order as u
% The file is overwritten every time.
t   = t(:);
pct = 100 * (u ./ u_max(:)).';                      % n x 3, order [ddx ddl ddth]
pct = pct(:, [1 3 2]);                              % -> trolley, slew, hoist

C = cell(numel(t) + 1, 5);
C(1, :) = {'Continuous', 'time_s', 'trolley_pct', 'slew_pct', 'hoist_pct'};
C(2:end, 2:5) = num2cell([t, pct]);

if isfile(file), delete(file); end                  % so no stale rows from a longer old file survive
writecell(C, file);
fprintf('Saved %d input rows to %s (peak |u|: trolley %.0f%%, slew %.0f%%, hoist %.0f%%)\n', ...
        numel(t), file, max(abs(pct), [], 1));
end
