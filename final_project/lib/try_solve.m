function [sol, ok, info] = try_solve(sol_guess, params, X0, kdir, opts)
% One fixed-tf bvp4c solve with every failure check. Returns the guess unchanged on failure.
    odefun = @(t, Z) hamiltonian_ode(t, Z, params);
    bcfun  = @(ya, yb) bc_loop(ya, yb, X0);
    lastwarn('');                                   % clear old warnings
    try
        sol = bvp4c(odefun, bcfun, sol_guess, opts);
    catch err                                       % e.g. "singular Jacobian"
        sol = sol_guess;  ok = false;  info = err.message;  return
    end
    msg = lastwarn;                                 % non-empty = bvp4c complained
    Zf  = sol.y(:, end);
    end_err = max(abs([Zf(1:10) - X0; Zf(11) - 2*pi*kdir]));
    ok   = isempty(msg) && sol.stats.maxres < 1e-3 && end_err < 1e-3;
    info = sprintf('maxres=%.1e end_err=%.1e mesh=%d %s', ...
                   sol.stats.maxres, end_err, numel(sol.x), msg);
    if ~ok, sol = sol_guess; end
end
