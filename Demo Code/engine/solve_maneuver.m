function S = solve_maneuver(cfg, warm)
% Behind-the-scenes solver: gets from wherever we are to the PMP solution of cfg.
%
%   1. Start: warm = last saved solution (fast: a few small steps), or a cold start
%      from a guessed loop with the pin pulled in where the trolley can reach it.
%   2. Geometry: continuation (homotopy) of start/finish/pin/limits from the old problem
%      to the new one, halving the step on any failed solve.
%   3. Time: fixed t_f -> continue t_f to the requested value.
%            'min'     -> free t_f with H(t_f) = 0, then raise the time price W_t until
%                         t_f stops dropping (min t_f + (1/W_t)∫L dt -> minimum time),
%                         then shrink a fixed t_f until a limit or the pin blocks it.
%   Any soft wall that ends up crossed is hardened (k -> 3k) and the step re-solved.
clock = tic;
if nargin > 1 && usable(warm, cfg)
    fprintf('Warm start from the solution saved %s (t_f = %.2f s)\n', warm.created, warm.t_f);
    sol = remesh(warm.sol, cfg.warm_mesh);  cur = adopt(warm.cfg, cfg);
else
    [sol, cur] = cold_start(cfg);
end

% --- 2. Geometry and limits, in the current time mode ---
geo = {'X_start', 'X_end', 'pin_x', 'pin_y', 'r_p', 'c_clr', 'z_top', 'l_off', 'sweep', ...
       'dxmax', 'dlmax', 'dthmax', 'ddxmax', 'ddlmax', 'ddthmax'};
to = cur;
for f = geo, to.(f{1}) = cfg.(f{1}); end
[sol, cur] = continuation(sol, cur, to, geo, 'geometry');

% --- 3. Time ---
if ~cfg.free_tf
    if cur.free_tf                                    % same problem, t_f just stops being unknown
        cur.t_f = sol.parameters;  cur.free_tf = false;  cur.W_t = 0;  sol = rmfield(sol, 'parameters');
    end
    to = cur;  to.t_f = cfg.t_f;
    [sol, cur] = continuation(sol, cur, to, {'t_f'}, 't_f');
else
    if ~cur.free_tf                                   % dJ*/dt_f = H(t_f): this W_t makes the current t_f optimal
        c = cur;  c.free_tf = true;
        c.W_t = max(-hamiltonian_value(sol.y(:, end), crane_params(cur)), 1e-4);
        [sol, ok, info] = bvp_attempt(sol, c);
        if ~ok, error('Switching to free t_f failed: %s', info); end
        cur = c;  cur.t_f = sol.parameters;
        fprintf('  free t_f      W_t = %8.3g  %s\n', cur.W_t, info);
    end
    [sol, cur] = min_time(sol, cur);
    [sol, cur] = squeeze_tf(sol, cur);
end

S = package(sol, cur);
fprintf('Solved in %.0f s: t_f = %.2f s, %s\n', toc(clock), S.t_f, ternary(S.M.admissible, 'all limits OK', ['VIOLATES: ' strjoin(S.M.why, ', ')]));
end

%% ---------------------------------------------------------------------------
function ok = usable(warm, cfg)
% A warm start only helps if it solved the same kind of problem on the same side of the pin.
ok = isstruct(warm) && isfield(warm, 'cfg') && all(isfield(warm.cfg, param_names())) ...
     && isfield(warm.cfg, 'sweep') && abs(warm.cfg.sweep - cfg.sweep) < pi;
if ~ok && isstruct(warm), fprintf('Saved solution is a different maneuver -> cold start\n'); end
end

function s = remesh(sol, n)
% Same solution on a coarse uniform mesh. The saved mesh is sized for the OLD problem's
% fine features; starting coarse lets bvp4c refine only where the new one needs it.
if isfield(sol, 'parameters')
    s = bvpinit(linspace(0, 1, n), @(tau) deval(sol, tau), sol.parameters);
else
    s = bvpinit(linspace(0, 1, n), @(tau) deval(sol, tau));
end
end

function c = adopt(old, new)
% Old problem with the new settings' bookkeeping (task, pin slew, loops, hardware) but the
% old geometry and hardened weights, which continuation moves to the new values.
c = old;
for f = {'task', 'tune', 'pin_slew_deg', 'passes', 'direction', 'hw', 'bvp', 'maxres', 'step0', 'step_min', ...
         'max_hard', 'hard_fac', 'W_t_max', 'W_t_fac0', 'tf_tol', 'tol', 'warm_mesh', 'drive_tau'}
    c.(f{1}) = new.(f{1});
end
end

function [sol, cur] = cold_start(cfg)
% Solve an easy version first: fixed t_f from the guess and the pin pulled in until the
% guessed path is inside the trolley's reach. Then push the pin out to where it really is.
fprintf('Cold start\n');
easy = cfg;  easy.free_tf = false;  easy.W_t = 0;
far  = easy.pin_x + easy.r_p + easy.c_clr + easy.guess_margin;
if far > easy.xmax - 0.01
    easy.pin_x = easy.pin_x - (far - (easy.xmax - 0.01));
    easy.sweep = target_sweep(easy);
end
[guess, easy.t_f] = initial_guess(easy);
[sol, ok, info] = bvp_attempt(guess, easy);
fprintf('  first solve   pin at %.3f m  t_f = %.1f s  %s\n', easy.pin_x, easy.t_f, info);
if ~ok, error('Cold start failed. Try a smaller clearance, a pin closer to the tower, or a fixed t_f.'); end
to = easy;  to.pin_x = cfg.pin_x;  to.sweep = cfg.sweep;
[sol, cur] = continuation(sol, easy, to, {'pin_x', 'sweep'}, 'pin');
end

function [sol, cur] = continuation(sol, from, to, fields, label)
% Homotopy cfg(λ) = (1-λ)·from + λ·to over the listed fields, λ: 0 -> 1, adaptive steps.
cur = from;
if isequal(cellfun(@(f) from.(f), fields, 'UniformOutput', false), ...
           cellfun(@(f) to.(f),   fields, 'UniformOutput', false)), return; end
lam = 0;  step = from.step0;
while lam < 1
    l2 = min(1, lam + step);
    c  = cur;
    for f = fields, c.(f{1}) = (1 - l2)*from.(f{1}) + l2*to.(f{1}); end
    [s2, ok, info] = bvp_attempt(sol, c);
    if ok, [s2, c, info] = harden(s2, c, info, false); end
    if ok
        sol = s2;  cur = c;  lam = l2;  step = min(2*step, 1);
        fprintf('  %-12s %5.1f %%  %s\n', label, 100*lam, info);
    else
        step = step/2;
        if step < from.step_min
            error('Stuck at %.1f %% of the %s change (%s). Last good solution kept; try a smaller change.', 100*lam, label, info);
        end
    end
end
end

function [sol, c, info] = harden(sol, c, info, strict)
% Penalty method: if a soft wall was crossed, stiffen it and re-solve from this solution.
% strict = true: return ok only if the result respects every limit.
for n = 1:c.max_hard
    M = metrics_of(sol, c);
    if M.admissible || any(strcmp(M.why, 'wrong side of pin')), break; end
    c2 = c;
    if any(strcmp(M.why, 'velocity')),                               c2.k_v   = c.hard_fac*c.k_v;   end
    if any(strcmp(M.why, 'x')),                                      c2.k_x   = c.hard_fac*c.k_x;   end
    if any(strcmp(M.why, 'payload-pin') | strcmp(M.why, 'cable-pin')), c2.W_obs = c.hard_fac*c.W_obs; end
    [s2, ok, info2] = bvp_attempt(sol, c2);
    if ~ok, break; end
    sol = s2;  c = c2;  info = [info2 sprintf('  (hardened: %s)', strjoin(M.why, ', '))];
end
if strict && ~metrics_of(sol, c).admissible, info = ''; end
end

function [sol, cur] = min_time(sol, cur)
% Raise W_t (price per second) until t_f converges. Only admissible steps are accepted,
% since a high time price will happily trade the soft walls for speed.
fac = cur.W_t_fac0;  n_small = 0;  n_fail = 0;
while cur.W_t < cur.W_t_max
    c = cur;  c.W_t = min(cur.W_t*fac, cur.W_t_max);
    [s2, ok, info] = bvp_attempt(sol, c);
    if ok, [s2, c, info] = harden(s2, c, info, true); ok = ~isempty(info); end
    if ok
        gain = 1 - s2.parameters / sol.parameters;
        sol = s2;  cur = c;  cur.t_f = sol.parameters;
        fprintf('  min time     W_t = %8.3g  %s  (-%.2f %%)\n', cur.W_t, info, 100*gain);
        if gain < cur.tf_tol, n_small = n_small + 1; else, n_small = 0; end
        if n_small >= 2, break; end
        fac = min(1.5*fac, 8);
    else
        fac = sqrt(fac);  n_fail = n_fail + 1;
        if fac < 1.02 || n_fail >= 2, break; end     % stalled: squeeze_tf takes over
    end
end
end

function [sol, cur] = squeeze_tf(sol, cur)
% Finish the minimum-time search with t_f fixed: shrink t_f until no admissible solution
% exists. (The W_t route can stall on an ill-conditioned Newton matrix well before the
% limits bind; this one only stops when a wall or the pin does.)
cur.t_f = sol.parameters;  cur.free_tf = false;  cur.W_t = 0;  sol = rmfield(sol, 'parameters');
d = 0.04;
while d > cur.tf_tol
    c = cur;  c.t_f = cur.t_f*(1 - d);
    [s2, ok, info] = bvp_attempt(sol, c);
    if ok, [s2, c, info] = harden(s2, c, info, true); ok = ~isempty(info); end
    if ok
        sol = s2;  cur = c;  d = min(1.5*d, 0.08);
        fprintf('  shrink t_f   %s\n', info);
    else
        d = d/2;
    end
end
end

function M = metrics_of(sol, c)
M = crane_metrics(crane_trajectory(package_raw(sol, c), linspace(0, tf_of(sol, c), 600)), c);
end

function S = package_raw(sol, c)
S = struct('cfg', c, 'sol', sol, 't_f', tf_of(sol, c));
end

function S = package(sol, c)
S = package_raw(sol, c);
c.t_f = S.t_f;  S.cfg = c;
S.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
S.M = crane_metrics(crane_trajectory(S, linspace(0, S.t_f, 2000)), c);
end

function tf = tf_of(sol, c)
if c.free_tf, tf = sol.parameters; else, tf = c.t_f; end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
