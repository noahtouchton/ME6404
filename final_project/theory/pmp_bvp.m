function B = pmp_bvp(cfg)
% The two-point boundary value problem that Pontryagin's principle hands bvp4c.
%
% Unknowns: Z = [x; λ] (11 states + 11 costates) on τ = t/t_f ∈ [0, 1], plus t_f itself
% when the final time is free. EOM.m derived everything that goes in here:
%
%   ODEs   dZ/dτ = t_f · F(Z),   F = [ ∂H/∂λ ; -∂H/∂x ] with u = u*(x, λ) substituted
%
%   BCs    x(0)   = x_start     10   start at rest at the start pose
%          η(0)   = 0            1   winding is counted from the start
%          x(t_f) = x_end       10   end at rest at the finish pose
%          λ_η(t_f) = 0          1   η(t_f) is NOT fixed: it only depends on the path's
%                                     endpoints and topology (gauge), so fixing it would
%                                     make the Newton matrix singular. The loop count comes
%                                     from the initial guess + the pin wall, and is checked after.
%          H(t_f) = 0            1   only if t_f is free (transversality: J doesn't depend on t_f
%                                     directly, so the optimal t_f has dJ/dt_f = H(t_f) = 0)
%
% Returns function handles for bvp4c, with the (cheap, exact) Jacobian of the BCs.
% bvp4c finite-differences the ODE Jacobian itself, one vectorized call per column
% (measured much faster than evaluating the symbolic 22x22 Jacobian point by point).
p  = crane_params(cfg);
F  = @(Z) [hamiltonian_ode(Z, p);  zeros(1, size(Z, 2))];   % 22nd row: λ_η' = -∂H/∂η = 0 (η is cyclic)
                                                           % Z may hold every mesh point at once (vectorized)
X0 = cfg.X_start;
Xf = cfg.X_end;

bc  = @(ya, yb) [ya(1:10) - X0;  ya(11);  yb(1:10) - Xf;  yb(22)];
Ia  = [eye(11), zeros(11); zeros(11, 22)];            % ∂bc/∂ya
Ib  = [zeros(11, 22); eye(10), zeros(10, 12); zeros(1, 21), 1];   % ∂bc/∂yb

B.free_tf = cfg.free_tf;
if cfg.free_tf
    B.ode    = @(tau, Z, tf) tf * F(Z);
    B.bc     = @(ya, yb, tf) [bc(ya, yb);  hamiltonian_value(yb, p)];
    B.bcjac  = @(ya, yb, tf) deal([Ia; zeros(1, 22)], [Ib; hamiltonian_grad(yb, p)], zeros(23, 1));
else
    tf = cfg.t_f;
    B.ode    = @(tau, Z) tf * F(Z);
    B.bc     = bc;
    B.bcjac  = @(ya, yb) deal(Ia, Ib);
end
end
