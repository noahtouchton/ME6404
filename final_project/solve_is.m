%% solve_is.m - Input-shaped loops around the pin
% The same task as solve_oc.m, built instead from three point-to-point moves, each made of
% acceleration ramps at the hardware limits convolved with an input shaper (impulses whose
% swings cancel each other):
%   stage : slew to -th_s while the trolley runs out to xmax
%   sweep : slew from -th_s to +th_s; slewing flings the payload outward, past xmax, as it
%           crosses the pin line
%   return: trolley back to x0 and slew back to 0; the payload passes inboard of the pin
% Designs, each checked on the full nonlinear model and saved in results/:
%   classic : a standard shaper on every move (ZV, ZVD or EI). Shaper type, sweep speed and
%             sweep length come from a scan on a fast linear model. Shaping each axis alone
%             ignores how slewing couples them, so some swing is left.     -> is_classic.mat
%   cancel  : ZV sweep; the stage and return moves get 2-impulse shapers solved on the full
%             model so the payload is at rest when the sweep starts and at the end. Each of
%             those moves also tries a few slew speeds and trolley-slew delays. -> is_cancel.mat
%   robust  : like cancel, with a ZVD sweep and 3-impulse shapers that make the swing as small
%             as possible at cable lengths 0.95, 1.00 and 1.05 m at once.  -> is_robust.mat
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'model', 'generated'));

%% Settings
designs = {'classic', 'cancel', 'robust'};       % which designs to make
o.min_margin   = 0.010;                  % payload and cable at least 10 mm outside the pin clearance
o.shaper_types = {'ZV', 'ZVD', 'EI'};    % standard shapers the classic design tries
o.rate_grid    = 0.75:0.025:1;           % sweep slew speed, fraction of dthmax
o.coast_grid   = 0:0.05:2;               % full-speed part of the sweep, in swing periods
o.move_rates   = [1, 0.8, 0.6];          % stage and return slew speeds to try
o.stage_delay  = [0, 0.5, 1, 1.5];       % stage: the trolley starts this long after the slew (s)
o.return_delay = [0, 0.5, 1, 1.5, 2, 3]; % return: the slew starts this long after the trolley (s)

cfg = crane.task();
for name = designs
    t0 = tic;
    switch name{1}
        case 'classic', [d, label] = design_classic(cfg, o);
        case 'cancel',  [d, label] = design_cancel(cfg, o, 'ZV', 2, 0);
        case 'robust',  [d, label] = design_cancel(cfg, o, 'ZVD', 3, 0.05);
    end
    if isempty(d), continue; end
    [cmd, info] = build_loop(cfg, d);
    S = struct('cfg', cfg, 'kind', 'commands', 'cmd', cmd, 'design', d, 'info', info, 'label', label, ...
               't_f', cmd.t_end, 'source', ['solve_is ' name{1}], 'created', char(datetime('now')), ...
               'compute_time', toc(t0));
    crane.report(label, crane.metrics(crane.trajectory(S, linspace(0, S.t_f, 4000)), cfg), cfg);
    save(crane.file(['is_' name{1} '.mat']), '-struct', 'S');
end

%% Designs
function [d, label] = design_classic(cfg, o)
% Standard shaper on every move. Scan shaper type x sweep speed x sweep length on the linear
% model and keep the loop with the least swing left that clears the pin and loops once.
    w = sqrt(cfg.g / cfg.l0);
    best = struct('resid', inf);
    for k = 1:numel(o.shaper_types)
        sh = crane.shaper(o.shaper_types{k}, w);
        for rate = o.rate_grid
            for c = o.coast_grid
                dk = loop_design(cfg, sh, rate, c);
                g  = ground_model(build_loop(cfg, dk), cfg);
                ok = g.pin_margin >= o.min_margin && g.cable_margin >= o.min_margin && abs(g.windings - cfg.kdir) < 0.25;
                if ok && g.resid < best.resid
                    best = struct('resid', g.resid, 'd', dk, 'type', o.shaper_types{k}, 'rate', rate, 'coast', c);
                end
            end
        end
    end
    d = [];  label = '';
    if isinf(best.resid), fprintf('classic: no loop meets the pin margin\n'); return; end
    d = best.d;  label = sprintf('IS classic (%s)', best.type);
    fprintf('classic: %s shapers, sweep at %.1f %% speed, %.2f periods at full speed (linear model: %.1f mm left)\n', ...
        best.type, 100*best.rate, best.coast, 1000*best.resid);
end

function [d, label] = design_cancel(cfg, o, sweep_type, n_imp, dl)
% Standard shaper on the sweep. The stage and return shapers (n_imp impulses per axis) are
% solved on the full model so the payload is at rest after each, at cable lengths
% l0*(1 + [-dl 0 dl]) (dl = 0: exactly, at l0 only).
    w = sqrt(cfg.g / cfg.l0);  Tp = 2*pi/w;
    robust = dl > 0;
    l_set = cfg.l0*(1 + unique([-dl, 0, dl]));
    sh = crane.shaper(sweep_type, w);
    q0 = shaper_params_init(sh, n_imp);  q0 = [q0; q0];          % trolley and slew shapers start standard
    [lb, ub] = param_bounds(n_imp, Tp);  lb = [lb; lb];  ub = [ub; ub];
    lmo = struct('tol', 1e-8, 'max_iter', 30);
    if robust, label = sprintf('IS cancel robust (%s)', sweep_type); else, label = sprintf('IS cancel (%s)', sweep_type); end

    % 1. Sweep: the fastest speed and length that clear the pin, starting from rest (linear model)
    d = [];  best_t = inf;
    for rate = o.rate_grid
        for c = o.coast_grid(o.coast_grid <= 1)
            dc = loop_design(cfg, sh, rate, c);
            [cmd, info] = build_loop(cfg, dc);
            g = ground_model(sub_command(cmd, 2), cfg, [dc.x_top; cfg.l0; -info.th_s]);
            if g.pin_margin >= o.min_margin && g.cable_margin >= o.min_margin && cmd.t_end < best_t
                best_t = cmd.t_end;  d = dc;
            end
        end
    end
    if isempty(d), fprintf('%s: no sweep meets the pin margin\n', label); return; end

    % 2. Stage: for each slew speed and trolley delay, solve its shapers so the payload is at
    %    rest when it ends; keep the fastest that works
    R = [];  Q = {};
    for rs = o.move_rates
        for dly = o.stage_delay
            d.rate_stage = rs;  d.t_stage_trolley = dly;
            [q, F] = solve_lm(@(q) stage_resid(cfg, d, q, n_imp, l_set), q0, lb, ub, lmo);
            [~, info] = build_loop(cfg, set_shapers(d, 'stage', q, n_imp));
            R(end+1, :) = [rms_swing(F), info.stage(2), rs, dly];     %#ok<AGROW>
            Q{end+1} = q;                                              %#ok<AGROW>
        end
    end
    i = choose(R, robust);
    d.rate_stage = R(i,3);  d.t_stage_trolley = R(i,4);  d = set_shapers(d, 'stage', Q{i}, n_imp);

    % 3. Return: simulate up to its start (per cable length), then the same search
    [cmd, info] = build_loop(cfg, d);
    t_ret = info.sweep(2);
    X_ret = zeros(11, numel(l_set));
    for k = 1:numel(l_set)
        X0 = [cfg.X0; 0];  X0(2) = l_set(k);
        T = crane.simulate(crane.command(cmd), [0, t_ret], cfg, X0);
        X_ret(:, k) = T.X(:, end);
    end
    d.t_return_trolley = t_ret;
    R = [];  Q = {};
    for rr = o.move_rates
        for dly = o.return_delay
            d.rate_return = rr;  d.t_return_slew = t_ret + dly;
            [q, F] = solve_lm(@(q) return_resid(cfg, d, q, n_imp, X_ret, t_ret), q0, lb, ub, lmo);
            cmd = build_loop(cfg, set_shapers(d, 'return', q, n_imp));
            R(end+1, :) = [rms_swing(F), cmd.t_end, rr, dly];          %#ok<AGROW>
            Q{end+1} = q;                                              %#ok<AGROW>
        end
    end
    i = choose(R, robust);
    d.rate_return = R(i,3);  d.t_return_slew = t_ret + R(i,4);  d = set_shapers(d, 'return', Q{i}, n_imp);
    fprintf('%s: sweep +/-%.1f deg, stage slew %.0f %%, return slew %.0f %%, swing left %.4f mm (full model)\n', ...
        label, rad2deg(d.rate*cfg.dthmax*d.hold/2), 100*d.rate_stage, 100*d.rate_return, 1000*R(i,1));
end

function i = choose(R, robust)
% Grid point to keep; rows of R are [swing left, end time, slew speed, delay].
% exact: the fastest that leaves no swing (< 0.01 mm), else the least swing.
% robust: the fastest within 20 % of the least swing over the cable lengths.
    if robust
        k = find(R(:,1) <= 1.2*min(R(:,1)));
    else
        k = find(R(:,1) < 1e-5);
        if isempty(k), [~, i] = min(R(:,1)); return; end
    end
    [~, j] = min(R(k,2));
    i = k(j);
end

function r = stage_resid(cfg, d, q, n_imp, l_set)
% Swing left after the stage move (full model), stacked over the cable lengths
    cmd = sub_command(build_loop(cfg, set_shapers(d, 'stage', q, n_imp)), 1);
    c = crane.command(cmd);
    r = zeros(4*numel(l_set), 1);
    for k = 1:numel(l_set)
        X0 = [cfg.X0; 0];  X0(2) = l_set(k);
        T = crane.simulate(c, [0, cmd.t_end], cfg, X0);
        r(4*k-3:4*k) = swing_state(T.X(:, end), cfg.g);
    end
end

function r = return_resid(cfg, d, q, n_imp, X_ret, t_ret)
% Swing left at the end (full model), simulating only the return from its start state
    cmd = build_loop(cfg, set_shapers(d, 'return', q, n_imp));
    c = crane.command(cmd);
    r = zeros(4*size(X_ret, 2), 1);
    for k = 1:size(X_ret, 2)
        T = crane.simulate(c, [t_ret, cmd.t_end], cfg, X_ret(:, k));
        r(4*k-3:4*k) = swing_state(T.X(:, end), cfg.g);
    end
end

function r = swing_state(X, g)
% Swing as 4 lengths (m): radial and tangential offset, and rate/w
    l = X(2);  w = sqrt(g/l);
    r = l*[X(4); X(9)/w; X(5); X(10)/w];
end

function s = rms_swing(F)
% rms swing amplitude over the cable lengths of a stacked residual
    s = sqrt(sum(F.^2) / (numel(F)/4));
end

function d = set_shapers(d, move, q, n_imp)
% Trolley and slew shapers of a move ('stage' or 'return') from parameters q
    n = numel(q)/2;
    d.(['sh_' move '_trolley']) = shaper_params(q(1:n), n_imp);
    d.(['sh_' move '_slew'])    = shaper_params(q(n+1:end), n_imp);
end

function sh = shaper_params(q, n_imp)
% Shaper with positive amplitudes summing to 1 from box-bounded parameters, so a plain
% least-squares solver can tune it: n_imp = 2: q = [A, D] -> [A, 1-A] at [0, D];
% n_imp = 3: q = [s1, s2, D2, D3] -> [s1, (1-s1)*s2, (1-s1)*(1-s2)] at [0, D2, D3]
    if n_imp == 2, sh = [q(1), 1 - q(1); 0, q(2)];
    else,          sh = [q(1), (1 - q(1))*q(2), (1 - q(1))*(1 - q(2)); 0, q(3), q(4)];
    end
end

function q = shaper_params_init(sh, n_imp)
% Parameters q that reproduce standard shaper sh (ZV for 2 impulses, ZVD or EI for 3)
    if n_imp == 2, q = [sh(1,1); sh(2,2)];
    else,          q = [sh(1,1); sh(1,2)/(1 - sh(1,1)); sh(2,2); sh(2,3)];
    end
end

function [lb, ub] = param_bounds(n_imp, Tp)
% Box for shaper_params: amplitudes stay positive, impulses within 1.5-2 swing periods
    if n_imp == 2, lb = [0; 0.05*Tp];   ub = [1; 1.5*Tp];
    else,          lb = [0; 0; 0; 0];   ub = [1; 1; 2*Tp; 2*Tp];
    end
end

%% The loop as a ramp command
function [cmd, info] = build_loop(cfg, d)
% Ramp command of the loop from the stage, sweep and return moves. d holds the sweep speed
% (rate, fraction of dthmax) and length (hold: time between its two velocity steps), x_top,
% each move's start time (NaN: when the previous move ends), each axis's shaper and optional
% stage and return slew speeds. cmd.part labels each ramp's move (1 stage, 2 sweep, 3 return).
    IX = 1;  ITH = 3;
    W = d.rate * cfg.dthmax;
    th_s = W * d.hold / 2;                                % the sweep goes from -th_s to +th_s
    rs = 1;  if isfield(d, 'rate_stage'),  rs = d.rate_stage;  end
    rr = 1;  if isfield(d, 'rate_return'), rr = d.rate_return; end
    R1 = crane.shape_steps(ITH, crane.move_steps(d.t_stage_slew, -th_s, rs*cfg.dthmax, cfg.ddthmax), d.sh_stage_slew, cfg.ddthmax);
    R2 = crane.shape_steps(IX,  crane.move_steps(d.t_stage_trolley, d.x_top - cfg.x0, cfg.dxmax, cfg.ddxmax), d.sh_stage_trolley, cfg.ddxmax);
    stage = span([R1; R2]);
    t_sw = d.t_sweep;  if isnan(t_sw), t_sw = stage(2); end
    R3 = crane.shape_steps(ITH, [t_sw, W; t_sw + d.hold, -W], d.sh_sweep, cfg.ddthmax);
    sweep = span(R3);
    t_rx = d.t_return_trolley;  if isnan(t_rx), t_rx = sweep(2); end
    t_rs = d.t_return_slew;     if isnan(t_rs), t_rs = sweep(2); end
    R4 = crane.shape_steps(IX,  crane.move_steps(t_rx, cfg.x0 - d.x_top, cfg.dxmax, cfg.ddxmax), d.sh_return_trolley, cfg.ddxmax);
    R5 = crane.shape_steps(ITH, crane.move_steps(t_rs, -th_s, rr*cfg.dthmax, cfg.ddthmax), d.sh_return_slew, cfg.ddthmax);
    cmd.ramps = [R1; R2; R3; R4; R5];
    cmd.t_end = max(cmd.ramps(:,2) + cmd.ramps(:,3));
    cmd.part  = [ones(size(R1,1) + size(R2,1), 1); 2*ones(size(R3,1), 1); 3*ones(size(R4,1) + size(R5,1), 1)];
    info.th_s  = th_s;
    info.stage = stage;
    info.sweep = sweep;
    info.return = span([R4; R5]);
end

function d = loop_design(cfg, sh, rate, coast)
% Default loop: shaper sh on every move, moves back to back, and a sweep at rate*dthmax that
% holds full speed for coast swing periods
    w = sqrt(cfg.g / cfg.l0);
    d.rate  = rate;
    d.hold  = sh(2,end) + rate*cfg.dthmax/cfg.ddthmax + coast*2*pi/w;
    d.x_top = cfg.xmax;
    d.t_stage_slew = 0;  d.t_stage_trolley = 0;
    d.t_sweep = NaN;  d.t_return_trolley = NaN;  d.t_return_slew = NaN;
    [d.sh_stage_slew, d.sh_stage_trolley, d.sh_sweep, d.sh_return_trolley, d.sh_return_slew] = deal(sh);
end

function c = sub_command(cmd, parts)
% Some moves of a loop command (1 stage, 2 sweep, 3 return) as a command of their own
    keep = ismember(cmd.part, parts);
    c.ramps = cmd.ramps(keep, :);  c.part = cmd.part(keep);
    c.t_end = max(c.ramps(:,2) + c.ramps(:,3));
end

function s = span(R)
% [first ramp start, last ramp end] of ramp rows R
    s = [min(R(:,2)), max(R(:,2) + R(:,3))];
end

%% Fast linear model and least squares
function G = ground_model(cmd, cfg, P0)
% Fast linear prediction of the payload for the design scans. Horizontal ground-frame model
% with the cable length fixed at l0: z'' = w^2 (p - z), where z is the payload and
% p = x e^{i theta} the trolley as complex numbers. It drops only swing-squared terms, so the
% centrifugal and Coriolis effects of slewing are in. P0 = [x; l; theta] where the actuators
% start, at rest, with the payload hanging still (default: the task's start). Returns the pin
% and cable margins (m), the windings, and the swing left at the end as resid (m).
    if nargin < 3, P0 = [cfg.x0; cfg.l0; 0]; end
    N  = max(ceil(cmd.t_end/2e-3), 1);  h = cmd.t_end / N;
    tt = linspace(0, cmd.t_end, 2*N + 1);                     % RK4 grid with half steps
    P  = crane.profile(cmd.ramps, tt, P0);
    p  = P(1,:) .* exp(1i*P(3,:));                            % trolley, ground frame
    w2 = cfg.g / cfg.l0;  w = sqrt(w2);
    z  = complex(zeros(1, N + 1));  zk = p(1);  vk = 0;  z(1) = zk;
    for k = 1:N
        p0 = p(2*k-1);  pm = p(2*k);  p1 = p(2*k+1);
        a1 = w2*(p0 - zk);
        z2 = zk + h/2*vk;  v2 = vk + h/2*a1;  a2 = w2*(pm - z2);
        z3 = zk + h/2*v2;  v3 = vk + h/2*a2;  a3 = w2*(pm - z3);
        z4 = zk + h*v3;    v4 = vk + h*a3;    a4 = w2*(p1 - z4);
        zk = zk + h/6*(vk + 2*v2 + 2*v3 + v4);
        vk = vk + h/6*(a1 + 2*a2 + 2*a3 + a4);
        z(k+1) = zk;
    end
    pg = p(1:2:end);  pin = cfg.R_p;
    G.pin_margin = min(abs(z - pin)) - (cfg.r_p + cfg.c_clr);
    c  = pg + (abs(cfg.z_top)/cfg.l0)*(z - pg);               % cable point at pin-top height
    AB = z - c;
    s  = min(max(real((pin - c).*conj(AB)) ./ max(abs(AB).^2, eps), 0), 1);
    G.cable_margin = min(abs(c + s.*AB - pin)) - cfg.r_p;
    ang = unwrap(angle(z - pin));
    G.windings = (ang(end) - ang(1)) / (2*pi);
    dz = z(end) - pg(end);                                    % actuators are at rest at the end
    G.resid = norm([real(dz), real(vk)/w, imag(dz), imag(vk)/w]);
end

function [q, F] = solve_lm(fun, q0, lb, ub, opts)
% Levenberg-Marquardt for min ||fun(q)||^2 with bounds lb <= q <= ub (projected steps) and a
% forward-difference Jacobian. Base MATLAB only: the Optimization Toolbox isn't installed.
    tol = opts.tol;  max_iter = opts.max_iter;  h = 1e-6;
    q = min(max(q0(:), lb), ub);
    F = fun(q);  F = F(:);  f = F.'*F;
    lambda = 1e-3;
    warn = warning('off', 'MATLAB:nearlySingularMatrix');    % a parameter at a bound can have no effect
    cleanup = onCleanup(@() warning(warn));
    for it = 1:max_iter
        if sqrt(f) <= tol, break; end
        J = zeros(numel(F), numel(q));
        for k = 1:numel(q)
            hk = h;  if q(k) + hk > ub(k), hk = -hk; end      % stay inside the box
            qk = q;  qk(k) = qk(k) + hk;
            Fk = fun(qk);
            J(:, k) = (Fk(:) - F) / hk;
        end
        g = J.'*F;  H = J.'*J;
        D = diag(H);  D = max(D, 1e-6*max(D) + eps);           % scaled damping, never zero
        improved = false;
        while lambda < 1e12
            qn = min(max(q - (H + lambda*diag(D)) \ g, lb), ub);
            Fn = fun(qn);  Fn = Fn(:);  fn = Fn.'*Fn;
            if fn < f
                step = norm(qn - q);
                q = qn;  F = Fn;  f = fn;  lambda = max(lambda/5, 1e-12);  improved = true;
                break
            end
            lambda = lambda*5;
        end
        if ~improved || step < 1e-12, break; end
    end
end
