%% EOM.m - Pontryagin's minimum principle for the tower crane, start to finish
% Run this once after changing the model or the cost. It writes the numeric
% functions the solver uses into ../generated/ (never edit those by hand).
%
%   1. Dynamics          x' = f(x, u)                   (Newton/transport theorem)
%   2. Cost              J  = ∫ L(x, u) dt              (time + effort + walls + pin)
%   3. Hamiltonian       H  = L + λᵀ f
%   4. Optimal control   ∂H/∂u = 0   ->  u* = u_max·tanh(-σ·u_max/r)
%   5. Costates          λ' = -∂H/∂x
%   6. Boundary values   x(0) = x_start, x(t_f) = x_end, λ_η(t_f) = 0, and H(t_f) = 0 if t_f is free
%                        (imposed in pmp_bvp.m, the only other file with theory in it)
clc; clear;

%% 1. Dynamics: spherical pendulum hanging from a trolley on a slewing jib
% Generalized coordinates: trolley radius x, hoist cable l, slew θ, swing angles φ (radial) and ψ (tangential).
% Inputs are the three drive accelerations u = [ẍ; l̈; θ̈].
syms g m T real
syms x dx ddx real
syms l dl ddl real
syms th dth ddth real
syms phi dphi ddphi real
syms psi dpsi ddpsi real
syms l_off real    % lab correction: if the measured swing period doesn't match the cable length, the
                   % pendulum effectively hangs at l + l_off (0 = trust the cable length)

q   = [x; l; th; phi; psi];
dq  = [dx; dl; dth; dphi; dpsi];
ddq = [ddx; ddl; ddth; ddphi; ddpsi];
fulldiff = @(expr) jacobian(expr, [q; dq]) * [dq; ddq];        % total time derivative (chain rule)

ax = [1; 0; 0];  ay = [0; 1; 0];  az = [0; 0; 1];               % A: jib frame (rotates with θ)
by = ay;
bz = -sin(phi)*ax + cos(phi)*az;                                % B: after the radial swing φ
bx = cross(by, bz);
cx = bx;                                                        % C: after the tangential swing ψ (cable along -cz)
cy = -sin(phi)*sin(psi)*ax + cos(psi)*ay + cos(phi)*sin(psi)*az;
cz = simplify(cross(cx, cy));

w_A   = dth*az;                                                 % jib angular velocity
ddt_G = @(v) fulldiff(v) + cross(w_A, v);                       % transport theorem, A -> ground

L_eff = l + l_off;
r_Q   = x*ax - L_eff*cz;                                        % swinging mass, from the tower base
a_Q   = simplify(ddt_G(simplify(ddt_G(r_Q))));
F_Q   = -m*g*az + T*cz;                                         % gravity + cable tension

% Newton along the two directions the tension has no component in -> two swing equations
eom_phi = simplify(dot(F_Q - m*a_Q, cx) / m) == 0;
eom_psi = simplify(dot(F_Q - m*a_Q, cy) / m) == 0;
sol_dd  = solve([eom_phi, eom_psi], [ddphi, ddpsi]);
ddphi_expr = simplify(sol_dd.ddphi);
ddpsi_expr = simplify(sol_dd.ddpsi);

% Payload position (the swinging mass, what must miss the pin), ground frame
r_Px = x + L_eff*sin(phi)*cos(psi);                             % in the jib frame
r_Py = L_eff*sin(psi);
r_Pz = -L_eff*cos(phi)*cos(psi);
Qp   = [cos(th) -sin(th); sin(th) cos(th)] * [r_Px; r_Py];      % horizontal, ground frame
dQp  = jacobian(Qp, q) * dq;

% Winding angle η around the pin: η' = rate the payload circles the pin center.
syms pin_x pin_y real
delta = Qp - [pin_x; pin_y];
d2    = delta.'*delta;                                          % squared payload-pin distance
w_eta = (delta(1)*dQp(2) - delta(2)*dQp(1)) / d2;

syms eta real
X = [x; l; th; phi; psi; dx; dl; dth; dphi; dpsi; eta];         % 11 states
U = [ddx; ddl; ddth];                                           % 3 inputs
f = [dx; dl; dth; dphi; dpsi; ddx; ddl; ddth; ddphi_expr; ddpsi_expr; w_eta];

%% 2. Cost  L = W_t + L_effort + L_walls + L_pin
syms W_t real                                           % price of one second (0 when t_f is fixed)
syms r1 r2 r3 real                                      % effort weights: trolley, hoist, slew
syms ddxmax ddlmax ddthmax real
syms xmin xmax lmin lmax dxmax dlmax dthmax real
syms k_x k_v real                                       % wall stiffness (penalty method: k -> ∞ recovers the hard limit)
syms r_p c_clr W_obs real                               % pin radius, clearance, pin-wall weight

% Effort: k·[u·atanh(u/u_max) + (u_max/2)·ln(1 - u²/u_max²)], chosen because its
% ∂/∂u = k·atanh(u/u_max), so H_u = 0 solves EXACTLY to a tanh: |u*| < u_max automatically.
U_max = [ddxmax; ddlmax; ddthmax];
r_vec = [r1; r2; r3];
k_u   = r_vec ./ U_max;
L_effort = sum(k_u .* (U .* atanh(U ./ U_max) + (U_max/2) .* log(1 - (U ./ U_max).^2)));

% Position/velocity limits as steep polynomial walls (≈0 inside, explode at the limit)
wall  = @(v, lo, hi, n) ((v - (lo+hi)/2) / ((hi-lo)/2))^n;
L_walls = 0.5*( k_x*wall(x, xmin, xmax, 20) + k_x*wall(l, lmin, lmax, 20) ...
              + k_v*(dx/dxmax)^8 + k_v*(dl/dlmax)^8 + k_v*(dth/dthmax)^8 );

% Pin: W_obs at the clearance circle, ~d⁻¹² inside it
L_pin = W_obs * ((r_p + c_clr)^2 / d2)^6;

L_cost = W_t + L_effort + L_walls + L_pin;

%% 3. Hamiltonian
syms px pl pth pphi ppsi pdx pdl pdth pdphi pdpsi peta real
P = [px; pl; pth; pphi; ppsi; pdx; pdl; pdth; pdphi; pdpsi; peta];   % costates λ
H = L_cost + P.'*f;

%% 4. Optimal control: interior minimum, ∂H/∂u = 0 (cone condition with no active bound)
% f is affine in u, so ∂(λᵀf)/∂u = σ (the switching function) doesn't depend on u:
%     ∂H/∂u = k·atanh(u/u_max) + σ = 0   ->   u* = u_max·tanh(-σ·u_max/r)
sigma = jacobian(P.'*f, U).';
assert(isempty(intersect(symvar(sigma), U)), 'σ must not depend on u (f affine in u)')
z     = -sigma .* U_max ./ r_vec;
U_opt = U_max .* tanh(z);
H_uu  = hessian(L_effort, U);       % Legendre-Clebsch: k/(u_max(1-u²/u_max²)) > 0 on the diagonal
assert(isequal(simplify(H_uu - diag(diag(H_uu))), sym(zeros(3))), 'H_uu should be diagonal')

%% 5. Costates: λ' = -∂H/∂x (u held fixed, then u* substituted)
dP_dt = -jacobian(H, X).';
assert(isequal(simplify(dP_dt(end)), sym(0)), 'λ_η should be constant: η never appears in H')

F = subs([f; dP_dt], U, U_opt);              % the 22 ODEs of the boundary value problem

% Hamiltonian along the optimum (rewritten so atanh(u*/u_max) never sees |u*/u_max| -> 1)
L_effort_opt = sum(r_vec .* (z .* tanh(z) - log(cosh(z))));
H_opt = subs(L_cost - L_effort, U, U_opt) + L_effort_opt + P.'*subs(f, U, U_opt);

%% 6. Export the numeric functions
% Parameter order is written to param_names.m, so engine/crane_params.m reads the
% names from there and the two can't drift out of sync.
params_sym = [g; l_off; r1; r2; r3; xmin; xmax; lmin; lmax; dxmax; dlmax; dthmax; ...
              ddxmax; ddlmax; ddthmax; pin_x; pin_y; r_p; c_clr; W_obs; k_x; k_v; W_t];
Z = [X; P];

gen = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'generated');
if ~isfolder(gen), mkdir(gen); end
fprintf('Writing generated functions to %s ...\n', gen);
% F(22) = λ_η' = 0 is left out (pmp_bvp.m adds it): a constant row would stop the function
% from evaluating the whole mesh in one call, which is what makes bvp4c fast.
matlabFunction(F(1:21),           'File', fullfile(gen, 'hamiltonian_ode'),   'Vars', {Z, params_sym});
matlabFunction(U_opt,             'File', fullfile(gen, 'optimal_inputs'),    'Vars', {Z, params_sym});
matlabFunction(H_opt,             'File', fullfile(gen, 'hamiltonian_value'), 'Vars', {Z, params_sym});
matlabFunction(jacobian(H_opt,Z), 'File', fullfile(gen, 'hamiltonian_grad'),  'Vars', {Z, params_sym});
matlabFunction([Qp; r_Pz],        'File', fullfile(gen, 'payload_position'),  'Vars', {X, params_sym});
matlabFunction(f,                 'File', fullfile(gen, 'crane_dynamics'),    'Vars', {X, U, params_sym});   % x' = f(x, u), for replaying commands

names = arrayfun(@char, params_sym, 'UniformOutput', false);
fid = fopen(fullfile(gen, 'param_names.m'), 'w');
fprintf(fid, 'function n = param_names()\n%% Order of the parameter vector in every generated function (written by theory/EOM.m).\nn = {%s};\nend\n', ...
        strjoin(strcat('''', names, ''''), ', '));
fclose(fid);
addpath(gen); rehash
fprintf('Done: %d states + %d costates, %d parameters.\n', numel(X), numel(P), numel(params_sym));
