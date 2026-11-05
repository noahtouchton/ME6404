function Y = loop_guess(t, t_f, X0, pin, kdir)
% Initial guess: payload travels |kdir| circles around the pin, cable straight down.
    tau = t / t_f;
    s   = tau - sin(2*pi*tau)/(2*pi);          % 0 -> 1, zero slope at both ends (starts/stops at rest)
    sd  = (1 - cos(2*pi*tau)) / t_f;           % ds/dt
    rho = norm([X0(1); 0] - pin);              % loop radius = start's distance from the pin
    a   = 2*pi*abs(kdir)*s;                    % angle around the pin: |kdir| full turns
    Q   = pin + rho*[-cos(a); -sign(kdir)*sin(a)];            % at a = 0 this is exactly the start point
    Qd  = rho*2*pi*abs(kdir)*sd*[sin(a); -sign(kdir)*cos(a)]; % its time derivative

    x   = norm(Q);                             % trolley radius (directly above the payload)
    th  = atan2(Q(2), Q(1));                   % jib angle
    dx  = dot(Q, Qd) / x;                      % radial rate
    dth = (Q(1)*Qd(2) - Q(2)*Qd(1)) / x^2;     % angular rate about the tower
    eta = kdir*2*pi*s;                         % winding angle so far

    Y = [x; X0(2); th; 0; 0; dx; 0; dth; 0; 0; eta; zeros(11,1)];   % costates filled in by caller
end
