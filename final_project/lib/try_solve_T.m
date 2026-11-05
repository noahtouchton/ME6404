function [sol, ok, info] = try_solve_T(sol_guess, params, X0, kdir, opts)
% Free-tf solve on tau in [0,1]; the unknown parameter sol.parameters is t_f.
% The 23rd BC is H(tf) = 0 (free final time, no terminal cost). Returns the guess on failure.
    odefun = @(tau, Z, T) T * hamiltonian_ode(T*tau, Z, params);
    bcfun  = @(ya, yb, T) [bc_loop(ya, yb, X0); hamiltonian_value(yb, params)];
    lastwarn('');
    try
        sol = bvp4c(odefun, bcfun, sol_guess, opts);
    catch err
        sol = sol_guess;  ok = false;  info = err.message;  return
    end
    msg = lastwarn;
    Zf  = sol.y(:, end);
    end_err = max(abs([Zf(1:10) - X0; Zf(11) - 2*pi*kdir]));
    ok   = isempty(msg) && sol.stats.maxres < 1e-3 && end_err < 1e-3 && sol.parameters > 0;
    info = sprintf('t_f=%.3f s maxres=%.1e end_err=%.1e mesh=%d %s', ...
                   sol.parameters, sol.stats.maxres, end_err, numel(sol.x), msg);
    if ~ok, sol = sol_guess; end
end
