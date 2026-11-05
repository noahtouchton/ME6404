clc; clear;

syms g m T real

syms x dx ddx real
syms l dl ddl real
syms th dth ddth real
syms phi dphi ddphi real
syms psi dpsi ddpsi real

q = [x; l; th; phi; psi];
dq = [dx; dl; dth; dphi; dpsi];
ddq = [ddx; ddl; ddth; ddphi; ddpsi];

%Total time derivative d/dt using Chain Rule
fulldiff = @(expr) jacobian(expr, [q;dq]) * [dq; ddq];

ax = [1; 0; 0];
ay = [0; 1; 0];
az = [0; 0; 1];

by = ay;
bz = -sin(phi)*ax + cos(phi)*az;
bx = cross(by, bz);

cx = bx;
cy = -sin(phi)*sin(psi)*ax + cos(psi)*ay + cos(phi)*sin(psi)*az;
cz = simplify(cross(cx, cy));

w_A = dth*az;
w_C = dth*az - dphi*by - dpsi*cx;

%Transport Theorum Helper for any vector expressed in A to G
ddt_G = @(v) fulldiff(v) + cross(w_A, v);

r_Q = x*ax - l*cz;

v_Q = simplify(ddt_G(r_Q));
a_Q = simplify(ddt_G(v_Q));

F_Q = -m*g*az + T*cz;

% Divide out m so MATLAB doesn't tack on '| m == 0'
eom_phi = simplify(dot(F_Q - m*a_Q, cx) / m) == 0;
eom_psi = simplify(dot(F_Q - m*a_Q, cy) / m) == 0;

disp('--- Raw Projected EOMs ---');
disp(eom_phi);
disp(eom_psi);

sol = solve([eom_phi, eom_psi], [ddphi, ddpsi]);
ddphi_expr = simplify(sol.ddphi);
ddpsi_expr = simplify(sol.ddpsi);

disp('--- State-Space Acceleration Equations ---');
fprintf('ddphi = \n'); pretty(ddphi_expr);
fprintf('ddpsi = \n'); pretty(ddpsi_expr);

% Hamiltonian

syms t real
syms px pl pth pphi ppsi pdx pdl pdth pdphi pdpsi real
syms r1 r2 r3 real
syms xmax lmax dxmax dlmax dthmax real
syms ddxmax ddlmax ddthmax real 
syms eta peta real
syms pin_x pin_y r_p c_clr W_obs k_x W_t real
syms k_v real   % velocity-wall weight: must stay >> W_t or speed beats the velocity limits

r_Qx = x + l*sin(phi)*cos(psi);
r_Qy = l*sin(psi);
r_rad = sqrt(r_Qx^2 + r_Qy^2);
r_z = -l*cos(phi)*cos(psi);

Q = [cos(th) -sin(th); sin(th) cos(th)]* [r_Qx; r_Qy];
dQ = jacobian(Q, q) * dq;

delta = Q - [pin_x; pin_y];
d2 = delta(1)^2 + delta(2)^2;

w = (delta(1)*dQ(2) - delta(2)*dQ(1))/d2;


X = [x; l; th; phi; psi; dx; dl; dth; dphi; dpsi; eta];
P = [px; pl; pth; pphi; ppsi; pdx; pdl; pdth; pdphi; pdpsi; peta];
U = [ddx; ddl; ddth];

f = [dx; dl; dth; dphi; dpsi; ddx; ddl; ddth; ddphi_expr; ddpsi_expr; w];

%Payload position in frame x (r and z)



% Smooth Barrier Functions for Inequality Constraints:
%    0 < x < xmax,  0 < l < lmax,  |dx| < dxmax,  |dl| < dlmax,  |dth| < dthmax
%    Using a 6th-order normalized polynomial wall: ((val - mid)/half)^6
% This prevents on off solutions by penalizing boundary cases more

B_x   = k_x * ((x-xmax/2)/(xmax/2))^20; % steep wall at x=0 and x=xmax
B_l   = ((l - lmax/2)/(lmax/2))^8;
B_dx  = k_v * (dx / dxmax )^8;
B_dl  = k_v * (dl / dlmax)^8;
B_dth = k_v * (dth / dthmax)^8;
L_barrier = 0.5 * (B_x + B_l + B_dx + B_dl + B_dth);

L_obs = W_obs * ((r_p + c_clr)^2/d2)^6; % W_obs at pin radius, grows fast inside of it. 
%L_control = 0.5 * (r1*(ddx/ddxmax)^2 + r2*(ddl/ddlmax)^2 + r3*(ddth/ddthmax)^2);

U_max = [ddxmax; ddlmax; ddthmax];
r_vec = [r1; r2; r3];
k_u = r_vec ./ U_max;
L_control = sum(k_u .* (U .* atanh(U ./ U_max) ...
    + (U_max/2) .* log(1 - (U ./ U_max).^2)));  % ∫ k·atanh(u/umax) du

L_cost = W_t + L_control + L_barrier + L_obs;

assert(~has(L_cost, t), 'L must not depend on t');

H = L_cost + dot(P, f);

%dH/dU = 0 solve for U*
sigma = jacobian(dot(P, f), U).';           % switching functions σ = ∂(λᵀf)/∂u
assert(isempty(intersect(symvar(sigma), U)), 'sigma should not depend on U (f is affine in U)')
z     = -sigma .* U_max ./ r_vec;           % tanh argument
U_opt = U_max .* tanh(z);                   % exact solution of H_u = 0

disp('d2L/dU2 (should be diagonal and positive):');
disp(simplify(hessian(L_control, U)))



% Condition 2: Costate Equations dP/dt = -dH/dX
dP_dt = -jacobian(H, X).';

% Substitute U_opt into State and Costate ODEs to get pure 20-ODE system in (X, P)
dX_dt_opt = subs(f, U, U_opt);
dP_dt_opt = subs(dP_dt, U, U_opt);
disp(dP_dt(end))     % λ̇_η: should print 0, since η appears nowhere in H
% Export to a fast MATLAB function for bvp4c!
params_sym = [g; r1; r2; r3; xmax; lmax; dxmax; dlmax; dthmax; ddxmax; ddlmax; ddthmax; pin_x; pin_y; r_p; c_clr; W_obs; k_x; W_t; k_v];   % 20 entries: order must match lib/crane_params.m
ode_22 = [dX_dt_opt; dP_dt_opt];
L_control_opt = sum(r_vec .* (z .* tanh(z) - log(cosh(z))));   % φ(u*) written in terms of z
H_opt = subs(L_cost - L_control, U, U_opt) + L_control_opt + dot(P, subs(f, U, U_opt));
Qpos = [Q; r_z];

gen = fullfile(fileparts(mfilename('fullpath')), 'generated');   % always write next to this file, whatever the current folder
if ~isfolder(gen), mkdir(gen); end
matlabFunction(ode_22, 'File', fullfile(gen, 'hamiltonian_ode'),   'Vars', {t, [X; P], params_sym});
matlabFunction(U_opt,  'File', fullfile(gen, 'optimal_inputs'),    'Vars', {[X; P], params_sym});
matlabFunction(H_opt,  'File', fullfile(gen, 'hamiltonian_value'), 'Vars', {[X; P], params_sym});
matlabFunction(Qpos,   'File', fullfile(gen, 'payload_position'),  'Vars', {X, params_sym});
addpath(gen); rehash                                     % make MATLAB see the new files
disp(size(hamiltonian_ode(0, rand(22,1), rand(20,1))))   % should print 22 1