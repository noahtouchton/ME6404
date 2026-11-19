classdef crane
% crane - Everything the scripts share, as static functions (call them as crane.task() etc.)
%   cfg = crane.task()                    the maneuver: hardware limits, pin, start = end state
%   p   = crane.params(cfg)               parameter vector for model/generated/ (EOM.m's order)
%   sh  = crane.shaper(type, w)           standard input shaper [amplitudes; times]
%   S   = crane.move_steps(t0, D, v, a)   velocity steps of a point-to-point move
%   R   = crane.shape_steps(ax, S, sh, a) those steps convolved with a shaper, as ramp rows
%   [P, V, A] = crane.profile(R, t, P0)   exact axis motion under ramp rows R
%   c   = crane.command(S)                actuator command of a saved solution or ramp rows
%   T   = crane.simulate(c, t, cfg, X0)   full nonlinear model under a command, sampled at t
%   T   = crane.trajectory(S, t)          states, inputs and geometry of a saved solution
%   M   = crane.metrics(T, cfg)           limit usage, pin clearance, leftover swing, verdict
%   crane.report(name, M, cfg)            short printout of M
%   f   = crane.file(name)                path of a file in results/
% A ramp row is [axis, t0, dur, acc] (axis 1 trolley, 2 hoist, 3 slew): acceleration acc on
% [t0, t0 + dur); overlapping rows add up. A saved solution S (results/*.mat) holds cfg,
% kind ('optimal' or 'commands'), t_f, label, source, created, compute_time, and either the
% bvp4c solution (sol, time_scaled: on tau = t/t_f) or the ramp command (cmd.ramps).
    methods (Static)

        function cfg = task()
        % The maneuver both methods solve, so they are always compared on the same task
            cfg.g       = 9.81;
            cfg.xmax    = 0.925;      cfg.lmax   = 1.450;      % m
            cfg.dxmax   = 72.6e-3;    cfg.dlmax  = 30.0e-3;    % m/s
            cfg.dthmax  = deg2rad(32.4);                       % rad/s
            cfg.ddxmax  = 72.6e-3;    cfg.ddlmax = 30.0e-3;    % m/s^2
            cfg.ddthmax = deg2rad(32.4);                       % rad/s^2
            cfg.R_p   = 0.90;     % pin center from the tower (m): clearing it takes 0.95 m > xmax
            cfg.r_p   = 0.03;     % pin radius (m)
            cfg.c_clr = 0.02;     % clearance: payload size + margin (m)
            cfg.z_top = -0.75;    % pin top relative to the jib (m)
            cfg.x0 = 0.70;  cfg.l0 = 1.00;                     % start = end, hanging at rest
            cfg.X0 = [cfg.x0; cfg.l0; zeros(8, 1)];
            cfg.kdir = 1;         % loops around the pin: +1 once counterclockwise
        end

        function p = params(cfg)
        % Order of params_sym in EOM.m. Cost weights exist only for optimal control; shaped
        % solutions get NaN there, which the dynamics never read.
            w = @(name) weight(cfg, name);
            p = [cfg.g; w('r1'); w('r2'); w('r3'); cfg.xmax; cfg.lmax; cfg.dxmax; cfg.dlmax; ...
                 cfg.dthmax; cfg.ddxmax; cfg.ddlmax; cfg.ddthmax; cfg.R_p; 0; cfg.r_p; cfg.c_clr; ...
                 w('W_obs'); w('k_x'); w('W_t'); w('k_v')];
        end

        function sh = shaper(type, w, Vtol)
        % Undamped input shaper for a mode at w (rad/s) as [amplitudes; times], amplitudes
        % summing to 1. ZV: 2 impulses over half a period. ZVD: 3 over one period, robust to
        % frequency error. EI: 3 over one period, allowing Vtol (default 5 %) vibration at w
        % for a wider insensitive band. 'none': a single unit impulse.
            if nargin < 3, Vtol = 0.05; end
            T = 2*pi/w;
            switch upper(type)
                case 'NONE', sh = [1; 0];
                case 'ZV',   sh = [1/2, 1/2; 0, T/2];
                case 'ZVD',  sh = [1/4, 1/2, 1/4; 0, T/2, T];
                case 'EI',   sh = [(1 + Vtol)/4, (1 - Vtol)/2, (1 + Vtol)/4; 0, T/2, T];
                otherwise,   error('crane.shaper: unknown type %s', type);
            end
        end

        function S = move_steps(t0, D, vmax, amax)
        % Velocity steps [time, dv] of the fastest unshaped move of signed distance D starting
        % at t0: up to speed v, then back to 0 after |D|/v, so the distance is exactly D
        % whatever the ramps look like. v = vmax, or less if the move is too short to reach it.
            if D == 0, S = zeros(0, 2); return; end
            v = min(vmax, sqrt(abs(D)*amax));
            S = [t0, sign(D)*v; t0 + abs(D)/v, -sign(D)*v];
        end

        function R = shape_steps(axis, S, sh, amax)
        % Ramp rows for velocity steps S = [time, dv] convolved with shaper sh: one copy of each
        % step per impulse, delayed by its time and scaled by its amplitude (a ramp of the same
        % length at amplitude*amax). Positive amplitudes summing to 1 never ask for more than
        % the unshaped move, and the final speed and position are those of the unshaped steps.
            R = zeros(0, 4);
            for k = 1:size(S, 1)
                for i = 1:size(sh, 2)
                    if sh(1,i) ~= 0
                        R(end+1, :) = [axis, S(k,1) + sh(2,i), abs(S(k,2))/amax, sh(1,i)*sign(S(k,2))*amax];   %#ok<AGROW>
                    end
                end
            end
        end

        function [P, V, A] = profile(R, t, P0)
        % Exact position, velocity and acceleration of the axes [x; l; theta] under ramp rows
        % R at times t (3 x numel(t) each), starting at rest at P0 (3x1)
            t = t(:).';
            P = repmat(P0(:), 1, numel(t));  V = zeros(3, numel(t));  A = zeros(3, numel(t));
            for k = 1:size(R, 1)
                ax = R(k,1);  d = R(k,3);  acc = R(k,4);
                s  = t - R(k,2);  sc = min(max(s, 0), d);       % time into this ramp, clipped
                A(ax,:) = A(ax,:) + acc*(s >= 0 & s < d);
                V(ax,:) = V(ax,:) + acc*sc;
                P(ax,:) = P(ax,:) + acc*(sc.^2/2 + d*max(s - d, 0));
            end
        end

        function c = command(S)
        % Command for crane.simulate: c.U(t) (3 x numel(t)), its breakpoints, and whether it is
        % piecewise constant. S is a saved solution or anything with ramp rows in S.ramps.
        % Optimal-control solutions replay U*(t) from the bvp4c costates, and zero after t_f.
            if isfield(S, 'ramps') || strcmp(S.kind, 'commands')
                if isfield(S, 'ramps'), R = S.ramps; else, R = S.cmd.ramps; end
                c.ramps  = R;
                c.U      = @(t) ramp_accel(R, t);
                c.breaks = unique([R(:,2); R(:,2) + R(:,3)]);
                c.pwc    = true;
            else
                p = crane.params(S.cfg);
                c.U      = @(t) optimal_command(S, p, t);
                c.breaks = S.t_f;
                c.pwc    = false;
            end
        end

        function T = simulate(c, t, cfg, X0)
        % Full nonlinear model (crane_dynamics) under command c, sampled at increasing times t,
        % from state X0 at t(1) (11x1, default: at rest at the start). Integrated piece by piece
        % between the command's breakpoints. Ramp commands also get their exact peak |accel|
        % and |velocity| per axis (U_peak, V_peak), which the samples could step over.
            if nargin < 4 || isempty(X0), X0 = [cfg.X0; 0]; end
            p = crane.params(cfg);
            opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
            t = t(:).';
            inner = c.breaks(c.breaks > t(1) & c.breaks < t(end));
            edges = unique([t(1); inner(:); t(end)]).';
            X = zeros(11, numel(t));  X(:, 1) = X0;  Xk = X0;
            for j = 1:numel(edges) - 1
                ta = edges(j);  tb = edges(j+1);
                if c.pwc
                    Uj  = c.U((ta + tb)/2);                      % constant on this piece
                    rhs = @(~, Xs) crane_dynamics(Xs, Uj, p);
                else
                    rhs = @(ts, Xs) crane_dynamics(Xs, c.U(ts), p);
                end
                sol = ode45(rhs, [ta tb], Xk, opts);
                in  = t >= ta & t <= tb;
                X(:, in) = deval(sol, t(in));
                Xk = sol.y(:, end);
            end
            T.t = t;  T.X = X;  T.U = c.U(t);  T.H = nan(1, numel(t));
            T.G = geometry(X, cfg);
            if c.pwc && numel(edges) > 1
                T.U_peak = max(abs(c.U((edges(1:end-1) + edges(2:end))/2)), [], 2);
                [~, V] = crane.profile(c.ramps, edges, zeros(3, 1));   % velocity peaks sit on edges
                T.V_peak = max(abs(V), [], 2);
            end
        end

        function T = trajectory(S, t)
        % Optimal control: read off the bvp4c solution, with U* and H from the costates.
        % Input shaping: only the command is stored, so the full model is simulated.
            if strcmp(S.kind, 'commands')
                T = crane.simulate(crane.command(S), t, S.cfg);
                return
            end
            p = crane.params(S.cfg);
            Z = solution_at(S, t);
            T.t = t;
            T.X = Z(1:11, :);
            T.U = optimal_inputs(Z, p);
            T.H = zeros(1, numel(t));
            for k = 1:numel(t), T.H(k) = hamiltonian_value(Z(:, k), p); end
            T.G = geometry(T.X, S.cfg);
        end

        function M = metrics(T, cfg)
        % M.admissible is true only if every hardware limit holds (0.1 % slack), the payload
        % stays outside the pin's clearance, the cable clears the pin and the payload goes
        % around it once; M.why lists what failed
            X = T.X;  G = T.G;
            v_pk = max(abs(X(6:8,:)), [], 2);
            u_pk = max(abs(T.U), [], 2);
            if isfield(T, 'V_peak'), v_pk = max(v_pk, T.V_peak); end
            if isfield(T, 'U_peak'), u_pk = max(u_pk, T.U_peak); end
            M.t_f         = T.t(end);
            M.H_spread    = max(T.H) - min(T.H);                % NaN for shaped solutions
            M.x_ratio     = max(X(1,:)) / cfg.xmax;
            M.v_ratio     = (v_pk ./ [cfg.dxmax; cfg.dlmax; cfg.dthmax]).';
            M.u_ratio     = (u_pk ./ [cfg.ddxmax; cfg.ddlmax; cfg.ddthmax]).';
            M.d_pin_min   = min(G.d_pin);
            M.d_cable_min = min(G.d_cable);
            M.reach_over  = max(G.r_rad) - cfg.xmax;
            M.windings    = X(11, end) / (2*pi);
            M.end_swing   = rad2deg(X(4:5, end)).';
            % Free swing left once the actuators stop (small-angle pendulum): amplitude
            % l*sqrt(angle^2 + (rate/w)^2), radial and tangential
            l = X(2, end);  w = sqrt(cfg.g / l);
            M.residual_parts = l * [hypot(X(4,end), X(9,end)/w), hypot(X(5,end), X(10,end)/w)];
            M.residual = hypot(M.residual_parts(1), M.residual_parts(2));
            tol = 1e-3;
            M.why = {};
            if M.x_ratio > 1 + tol,                           M.why{end+1} = 'x'; end
            if any(M.v_ratio > 1 + tol),                      M.why{end+1} = 'velocity'; end
            if any(M.u_ratio > 1 + tol),                      M.why{end+1} = 'acceleration'; end
            if M.d_pin_min < (cfg.r_p + cfg.c_clr)*(1 - tol), M.why{end+1} = 'payload-pin'; end
            if M.d_cable_min < cfg.r_p,                       M.why{end+1} = 'cable-pin'; end
            if abs(M.windings - cfg.kdir) > 0.25,             M.why{end+1} = 'windings'; end
            M.admissible = isempty(M.why);
        end

        function report(name, M, cfg)
        % Short printout of crane.metrics output (margins are beyond the required clearance)
            if M.admissible, verdict = 'ADMISSIBLE'; else, verdict = ['NOT ADMISSIBLE (' strjoin(M.why, ', ') ')']; end
            fprintf('%s: t_f = %.2f s | %s\n', name, M.t_f, verdict);
            fprintf('  reach past xmax %+6.1f mm | payload-pin margin %6.1f mm | cable-pin margin %6.1f mm\n', ...
                1000*M.reach_over, 1000*(M.d_pin_min - cfg.r_p - cfg.c_clr), 1000*(M.d_cable_min - cfg.r_p));
            fprintf('  |v|/vmax %.2f %.2f %.2f | |a|/amax %.2f %.2f %.2f | loops %+.3f\n', M.v_ratio, M.u_ratio, M.windings);
            fprintf('  swing left at t_f %.1f mm (radial %.1f, tangential %.1f)\n', 1000*M.residual, 1000*M.residual_parts);
        end

        function f = file(name)
        % Path of a file in results/ (created if missing)
            d = fullfile(fileparts(mfilename('fullpath')), 'results');
            if ~isfolder(d), mkdir(d); end
            f = fullfile(d, name);
        end
    end
end

%% Local functions
function G = geometry(X, cfg)
% Ground-frame positions (heights relative to the jib) and distances to the pin, for X (11 x N)
    p  = crane.params(cfg);
    N  = size(X, 2);  x = X(1, :);  th = X(3, :);
    G.Qp    = payload_position(X, p);                          % payload [X; Y; Z]
    G.P     = [x.*cos(th); x.*sin(th); zeros(1, N)];           % trolley
    G.Jib   = [cfg.xmax*cos(th); cfg.xmax*sin(th); zeros(1, N)];   % jib tip
    G.r_rad = hypot(G.Qp(1,:), G.Qp(2,:));                     % payload distance from the tower
    G.d_pin = hypot(G.Qp(1,:) - cfg.R_p, G.Qp(2,:));           % payload distance to the pin center
    % Cable vs pin: the part of the cable below the pin top must stay outside the pin
    G.d_cable = inf(1, N);
    low   = G.Qp(3,:) < cfg.z_top;                              % payload hangs below the pin top
    s_top = cfg.z_top ./ G.Qp(3, low);                          % cable fraction above pin-top height
    A     = G.P(1:2, low) + s_top .* (G.Qp(1:2, low) - G.P(1:2, low));   % cable point at pin-top height
    AB    = G.Qp(1:2, low) - A;
    C     = [cfg.R_p; 0];
    s     = max(0, min(1, sum((C - A).*AB, 1) ./ max(sum(AB.^2, 1), eps)));
    G.d_cable(low) = vecnorm(A + s.*AB - C);
end

function U = ramp_accel(R, t)
% Summed acceleration of ramp rows R at times t (3 x numel(t))
    t = t(:).';  U = zeros(3, numel(t));
    for k = 1:size(R, 1)
        on = t >= R(k,2) & t < R(k,2) + R(k,3);
        U(R(k,1), on) = U(R(k,1), on) + R(k,4);
    end
end

function U = optimal_command(S, p, t)
% U*(t) of an optimal-control solution from its costates, zero after t_f
    t = t(:).';  U = zeros(3, numel(t));
    k = t <= S.t_f;
    if any(k), U(:, k) = optimal_inputs(solution_at(S, t(k)), p); end
end

function Z = solution_at(S, t)
% bvp4c solution at real times t (free-t_f solutions live on tau = t/t_f)
    if S.time_scaled, Z = deval(S.sol, t / S.t_f); else, Z = deval(S.sol, t); end
end

function v = weight(cfg, name)
    if isfield(cfg, name), v = cfg.(name); else, v = NaN; end
end
