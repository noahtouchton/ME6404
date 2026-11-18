function [sol, ok, info] = bvp_attempt(guess, cfg)
% One bvp4c solve of the PMP boundary value problem for cfg, with every acceptance check.
% Returns the guess unchanged when the solve fails, so callers can just retry.
B = pmp_bvp(cfg);
opts = bvpset(cfg.bvp, 'BCJacobian', B.bcjac, 'Vectorized', 'on');

g = guess;                                   % match the guess to fixed / free t_f
if B.free_tf && ~isfield(g, 'parameters'), g.parameters = cfg.t_f; end
if ~B.free_tf && isfield(g, 'parameters'), g = rmfield(g, 'parameters'); end

lastwarn('');
w = warning('off', 'MATLAB:bvp4c:RelTolNotMet');      % checked below through maxres instead
try
    sol = bvp4c(B.ode, B.bc, g, opts);
    msg = lastwarn;
catch err
    warning(w);
    sol = guess;  ok = false;  info = strtok(err.message, newline);  return
end
warning(w);

if B.free_tf, tf = sol.parameters; else, tf = cfg.t_f; end
Zf = sol.y(:, end);
end_err   = max(abs(Zf(1:10) - cfg.X_end));
sweep_err = abs(Zf(11) - cfg.sweep);
ok = isempty(msg) && sol.stats.maxres < cfg.maxres && end_err < 1e-3 ...
     && sweep_err < 0.5 && tf > 0 && all(isfinite(sol.y(:)));
info = sprintf('t_f=%6.3f s  maxres=%.1e  mesh=%4d', tf, sol.stats.maxres, numel(sol.x));
if sweep_err >= 0.5, info = [info sprintf('  wrong side of pin (sweep %.0f° vs %.0f°)', rad2deg(Zf(11)), rad2deg(cfg.sweep))]; end
if ~isempty(msg),    info = [info '  ' strtok(msg, newline)]; end
if ~ok, sol = guess; end
end
