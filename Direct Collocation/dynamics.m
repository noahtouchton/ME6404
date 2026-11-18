function x_dot = dynamics(x, u)
    g = 9.81;
    % x = [x; l; th; phi; psi; dx; dl; dth; dphi; dpsi]
    % u = [ddx; ddl; ddth]

    [A, B] = crane_AB(x, g);

    x_dot = A + B*u;
end