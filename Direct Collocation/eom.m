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

X = [q; dq];
U = [ddx; ddl; ddth];
f = [dq; U; ddphi_expr; ddpsi_expr]; %f = x_dot so this comes up with the derivatives for each term

B = simplify(jacobian(f, U)); % taking df/dU will just give all the terms that are multiplied by U
A = simplify(subs(f, U, zeros(3,1))); %plugging in zero for U gives the A(x) terms

assert(isequal(simplify(f-(A+B*U)), sym(zeros(10,1))), 'f is not affine in u');

matlabFunction(A, B, 'File', 'crane_AB', 'Vars', {X, g}, 'Outputs', {'A', 'B'});

