%% lab_tune.m - Hot-tune the crane maneuver in the lab
% Edit sections 1-2, press Run. It warm-starts from the last solution (seconds, not
% minutes), checks the result, writes the continuous-velocity upload .xlsx into uploads/,
% replays that exact file through the crane model and plots a preview.
%
% The math lives in theory/ (EOM.m derives the minimum principle, pmp_bvp.m is the
% boundary value problem). engine/ does the numerical chores you shouldn't need to touch.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'theory'), fullfile(here, 'engine'), fullfile(here, 'export'), fullfile(here, 'generated'));

%% 1. The maneuver (what we care about)
% Poses are [trolley radius (mm), slew angle (deg), cable length (mm)], at rest.
task.start  = [700 0 1000];
task.finish = [700 0 1000];             % same as the start: loop around the pin and settle where it began

task.pin.radius_mm            = 30;     % pin's own radius
task.pin.radius_from_tower_mm = 900;    % where it stands: radius ...
task.pin.slew_deg             = 0;      % ... and slew angle of its center
task.pin.top_depth_mm         = 750;    % pin top, measured down from the jib (the cable can't cross the pin below this)

task.passes_behind = 1;                 % times the payload goes behind the pin (far side from the tower)
task.direction     = 'auto';            % 'auto' = shortest way, or '+' / '-' (sense of the slew angle)

task.t_f = 8;                           % seconds, or 'min' for the fastest admissible maneuver (~6.9 s here).
                                        % 'min' takes several minutes: find it before the lab, then use
                                        % that time + a margin here

%% 2. Lab corrections (change these based on what you see on the crane)
%  Symptom on the crane                                  Knob
%  ---------------------------------------------------   -------------------------------------------
%  Payload keeps swinging after the stop, and the        L_offset_mm  (or set measured_period_s and
%    swing period differs from the preview                            use the value it suggests)
%  Payload comes closer to the pin than the preview      clearance_mm up
%  Drives lag or saturate, path shape is off             limit_scale  down (e.g. 0.9)
%  An axis travels more/less than commanded              axis_gain    (export only: no re-solve)
tune.L_offset_mm  = 0;                  % pendulum effectively hangs at cable length + this
tune.clearance_mm = 20;                 % gap kept between payload and pin surface
tune.limit_scale  = 1.0;                % plan at this fraction of the velocity/acceleration limits
tune.axis_gain    = [1 1 1];            % [trolley slew hoist] actual / commanded travel

measured_period_s = [];                 % optional: residual swing period seen after the stop (s)

%% 3. Solve (skipped if sections 1-2 didn't change anything the optimizer sees)
sol_file = fullfile(here, 'results', 'solution.mat');
cfg  = build_problem(task, tune);
warm = [];
if isfile(sol_file), warm = load(sol_file); end
if same_problem(warm, cfg)
    S = warm;  fprintf('Problem unchanged: reusing the saved solution (t_f = %.2f s)\n', S.t_f);
else
    S = solve_maneuver(cfg, warm);
    if ~isfolder(fileparts(sol_file)), mkdir(fileparts(sol_file)); end
    save(sol_file, '-struct', 'S');
end
print_summary(S);

%% 4. Export the upload file and replay it through the model
out = fullfile(here, 'uploads', sprintf('trajectory_%s.xlsx', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
if ~isfolder(fileparts(out)), mkdir(fileparts(out)); end
cmd = export_continuous(S, out, tune);
save(strrep(out, '.xlsx', '.mat'), '-struct', 'S');       % the solution that made this file
R = simulate_commands(S, cmd);
fprintf('\nUpload file: %s\n', out);
fprintf('  replay of the file: %.2f s, end error %+.1f mm / %+.2f deg / %+.1f mm, pin gap %.1f mm, cable gap %.1f mm\n', ...
        R.t_cmd, R.end_pose_err, R.d_pin_min_mm - 1000*S.cfg.r_p, R.d_cable_min_mm - 1000*S.cfg.r_p);
fprintf('  residual swing after the stop: %.2f mm peak   (peak command %.0f / %.0f / %.0f %%)\n', ...
        R.residual_peak_mm, cmd.peak_pct);

%% 5. Suggestion from a measured swing period
if ~isempty(measured_period_s)
    L_eff = S.cfg.g * (measured_period_s/(2*pi))^2;
    fprintf('\nMeasured period %.2f s -> pendulum length %.0f mm -> set tune.L_offset_mm = %.0f\n', ...
            measured_period_s, 1000*L_eff, 1000*L_eff - task.finish(3));
end

%% 6. Preview
plot_preview(S, R);

%% ---------------------------------------------------------------------------
function same = same_problem(warm, cfg)
% True if the saved solution was solved for exactly these optimizer-visible settings.
same = false;
if ~isstruct(warm) || ~isfield(warm, 'cfg') || ~isfield(warm.cfg, 'task'), return; end
strip = @(t) rmfield(t, intersect(fieldnames(t), {'axis_gain'}));   % axis_gain only affects the export
same = isequal(warm.cfg.task, cfg.task) && isequal(strip(warm.cfg.tune), strip(cfg.tune)) ...
       && all(isfield(warm.cfg, param_names()));
end

function print_summary(S)
M = S.M;  c = S.cfg;
fprintf('\n  t_f %.2f s | swept %.0f deg around the pin | reach past trolley limit %+.1f mm\n', ...
        S.t_f, M.sweep_deg, 1000*M.reach_over);
fprintf('  limits used: v = %3.0f %3.0f %3.0f %%   a = %3.0f %3.0f %3.0f %%  (trolley hoist slew)\n', 100*M.v_ratio, 100*M.u_ratio);
fprintf('  payload-pin gap %.1f mm (asked %.0f) | cable-pin gap %.1f mm | max swing %.1f / %.1f deg\n', ...
        1000*(M.d_pin_min - c.r_p), 1000*c.c_clr, 1000*(M.d_cable_min - c.r_p), M.swing_max_deg);
fprintf('  PMP check: H spread %.1e (should be ~0)   %s\n', M.H_spread, ...
        ternary(M.admissible, 'ADMISSIBLE', ['VIOLATES: ' strjoin(M.why, ', ')]));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
