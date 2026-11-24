function dist = obstacle_sdf(p, ob, R_pay)
% Signed distance from a payload of radius R_pay centered at p (3x1, ground frame,
% Z up from the floor) to obstacle ob (one element of the obs struct from optimize.m).
% Positive = clear, negative = penetrating. Shared by optimize.m and sim_trajectory.m.
switch ob.type
    case 'box'
        % Vector from box center to payload center, in the box's frame (optional ob.yaw =
        % rotation about Z, rad), folded into the positive octant
        r = p - ob.center;
        if isfield(ob, 'yaw') && ~isempty(ob.yaw) && ob.yaw ~= 0
            c = cos(ob.yaw);  s = sin(ob.yaw);
            r = [c*r(1) + s*r(2); -s*r(1) + c*r(2); r(3)];
        end
        q = abs(r) - ob.half;
    case 'cylinder'
        d_radial   = norm(p(1:2) - ob.xy) - ob.R;
        d_vertical = abs(p(3) - 0.5*(ob.z(1) + ob.z(2))) - 0.5*(ob.z(2) - ob.z(1));
        q = [d_radial; d_vertical];
    otherwise
        error('obstacle_sdf: unknown obstacle type ''%s''', ob.type);
end
dist = (norm(max(q, 0)) + min(max(q), 0)) - R_pay;
end
