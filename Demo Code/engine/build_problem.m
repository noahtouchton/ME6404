function cfg = build_problem(task, tune)
% Lab units (mm, deg) -> model problem (SI, pin-aligned frame).
%
% Model frame: θ = 0 points at the pin, so the pin sits at (R_pin, 0) and
% θ_model = slew - pin slew. Every number the solver uses ends up in cfg.
H   = hidden_settings();
cfg = H;

% --- What we care about ---
cfg.pin_slew_deg = task.pin.slew_deg;
cfg.pin_x  = task.pin.radius_from_tower_mm / 1000;
cfg.pin_y  = 0;
cfg.r_p    = task.pin.radius_mm / 1000;
cfg.z_top  = -task.pin.top_depth_mm / 1000;       % pin top, relative to the trolley rail (down is negative)
cfg.passes    = task.passes_behind;
cfg.direction = task.direction;

pose = @(p) [p(1)/1000; p(3)/1000; deg2rad(p(2) - task.pin.slew_deg); zeros(7, 1)];  % [x; l; θ; φ ψ; rates] at rest
cfg.X_start = pose(task.start);
cfg.X_end   = pose(task.finish);

if ischar(task.t_f) || isstring(task.t_f)
    assert(strcmpi(task.t_f, 'min'), 't_f must be a number of seconds or ''min''');
    cfg.free_tf = true;   cfg.t_f = NaN;              % found by the solver
else
    cfg.free_tf = false;  cfg.t_f = task.t_f;  cfg.W_t = 0;
end

% --- Lab corrections ---
cfg.l_off = tune.L_offset_mm / 1000;
cfg.c_clr = tune.clearance_mm / 1000;
s = tune.limit_scale;                                 % plan inside the limits, leave margin for the real drives
cfg.dxmax  = s*H.dxmax;   cfg.dlmax  = s*H.dlmax;   cfg.dthmax  = s*H.dthmax;
cfg.ddxmax = s*H.ddxmax;  cfg.ddlmax = s*H.ddlmax;  cfg.ddthmax = s*H.ddthmax;
cfg.hw = H;                                           % unscaled hardware limits, for the final checks/export

% --- Target winding: total angle the payload sweeps around the pin ---
cfg.sweep = target_sweep(cfg);
if ~isfield(cfg, 'W_t'), cfg.W_t = 0; end
cfg.task = task;  cfg.tune = tune;
end
