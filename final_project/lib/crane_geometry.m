function G = crane_geometry(X, cfg)
% Ground-frame positions and pin distances for states X (11 x N).
% All heights are relative to the jib (z = 0 at the jib, payload z < 0).
params = crane_params(cfg);
N  = size(X, 2);
x  = X(1, :);
th = X(3, :);

G.Qp = zeros(3, N);                                % payload [X; Y; Z]
for k = 1:N
    G.Qp(:, k) = payload_position(X(:, k), params);
end
G.P   = [x.*cos(th); x.*sin(th); zeros(1, N)];                 % trolley
G.Jib = [cfg.xmax*cos(th); cfg.xmax*sin(th); zeros(1, N)];     % jib tip

G.r_rad = hypot(G.Qp(1,:), G.Qp(2,:));             % payload radius from the tower
G.d_pin = hypot(G.Qp(1,:) - cfg.R_p, G.Qp(2,:));   % payload distance to the pin center

% Cable vs pin: the part of the cable below the pin top must stay outside the pin
G.d_cable = inf(1, N);
for k = 1:N
    if G.Qp(3,k) < cfg.z_top                                    % payload hangs below the pin top
        s_top = cfg.z_top / G.Qp(3,k);                          % cable fraction above pin-top height
        A = G.P(1:2,k) + s_top*(G.Qp(1:2,k) - G.P(1:2,k));      % cable point at pin-top height
        G.d_cable(k) = seg_point_dist(A, G.Qp(1:2,k), [cfg.R_p; 0]);
    end
end
end

function d = seg_point_dist(A, B, C)
% Shortest distance from point C to the line segment A-B (all 2x1).
    AB = B - A;
    s  = max(0, min(1, dot(C - A, AB) / max(dot(AB, AB), eps)));   % closest point's fraction along AB
    d  = norm(A + s*AB - C);
end
