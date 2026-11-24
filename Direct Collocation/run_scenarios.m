%% run_scenarios.m - Five test maneuvers for comparing against input shaping
%
% HOW TO USE
%   1. Pick a scenario number in STEP 1.
%   2. Change the knobs for THAT scenario only, in its section in STEP 2. Every knob has
%      a comment saying what it does and what units it is in.
%   3. Press Run. The solver (Pontryagin / bvp4c, the same method as optimizeBVP.m) prints
%      its progress, then the 3D sim plays the result.
%   Not sure your obstacle or path setup makes sense? Set PREVIEW_ONLY = true first. It
%   draws the starting guess and the obstacles from above in about a second, without solving.
%
% UNITS AND COORDINATES
%   Lengths are in meters, angles in degrees.
%   The crane mast is at X = 0, Y = 0. Slew angle 0 deg points along +X, 90 deg along +Y.
%   A crane POSE is written [trolley, cable, slew]:
%     trolley = distance of the trolley from the mast   (0.02 to 0.905 m)
%     cable   = cable length below the jib               (0.2 to 1.43 m)
%     slew    = jib angle                                (any value, deg)
%   The jib is 2.0 m above the floor, so the payload hangs at a height of 2.0 - cable.
%   Obstacle positions are floor X/Y points in meters. To find the X/Y of a crane pose:
%   X = trolley*cosd(slew), Y = trolley*sind(slew).
%
% THE STARTING GUESS
%   The solver needs a rough path to start from. Each scenario builds one for you from its
%   knobs, so you rarely need to touch it. If a solve fails, the knobs worth trying first are
%   tf_guess (the guessed move time) and guess_clearance (how wide the guess goes around
%   obstacles).
%
% OUTPUT (in the results/ folder)
%   <name>.mat    the solved trajectory (same format as optimizeBVP.m)
%   <name>.xlsx   the inputs as % of max acceleration, for the crane GUI upload
%   The last solve is also copied to trajectory.mat, which sim_trajectory.m plays.
%   The hardware limits are fixed and set in solve_scenario_bvp.m.

clc; clear; close all;

%% ======================= STEP 1: choose what to run =======================
SCENARIO     = 1;       % 1 angle sweep | 2 pin | 3 wall | 4 lift while swinging | 5 trolley + slew
PREVIEW_ONLY = false;   % true = only draw the starting guess and obstacles, don't solve
RUN_SIM      = true;    % true = play the result in sim_trajectory.m after solving
SPEED_WEIGHT = 3.0;     % speed vs. smoothness. 3 = near minimum time (fast, hits the limits).
                        % Lower (e.g. 0.5) = slower, gentler moves. Use the same value for
                        % every scenario you plan to compare.

%% ======================= STEP 2: scenario knobs ===========================

% ---------- Scenario 1: ANGLE SWEEP ----------------------------------------
% A pure slew, rest to rest, at a fixed trolley position and cable length. One solve per
% angle in sweep_angles, so you can plot move time vs. slew angle. Short slews are where
% an input shaper's added delay hurts most.
S1.trolley      = 0.6;            % m, trolley distance from the mast (held for the whole move)
S1.cable        = 1.0;            % m, cable length (held). Longer = slower swing period
S1.start_slew   = 0;              % deg, where every move starts
S1.sweep_angles = [15 30 60 90];  % deg, how far to slew. Each one is a separate solve.
                                  % Negative = slew the other way.
S1.tf_guess     = [];             % s, guessed move time. [] = estimated automatically

% ---------- Scenario 2: PIN ------------------------------------------------
% A vertical pin (a thin post) the payload must go around.
%   mode = 'pass': move from start to finish, passing the pin on one side.
%   mode = 'loop': circle all the way around the pin and come back to the start pose.
S2.mode            = 'pass';
S2.pin_xy          = [0.6 0.0];   % m, floor X/Y of the pin center
S2.pin_radius      = 0.03;        % m, pin radius (the payload's own 8 cm radius is added for you)
S2.pin_height      = 1.5;         % m, pin height, used for drawing and the clearance check only
                                  % (the solver keeps the payload out of the pin's whole footprint)
S2.start           = [0.6 1.0 -40];   % pose [trolley m, cable m, slew deg]
S2.finish          = [0.6 1.0  40];   % pose, 'pass' mode only ('loop' ends at start)
S2.pass_side       = 'outside';   % 'pass' mode: 'outside' = go around the far side of the pin
                                  % (away from the mast), 'inside' = the side nearer the mast
S2.loop_direction  = 'ccw';       % 'loop' mode: 'ccw' or 'cw', seen from above
S2.guess_clearance = 0.08;        % m, how far from the pin's keep-out zone the guess path
                                  % passes. Bigger = safer start, possibly a slower answer
S2.tf_guess        = [];          % s, guessed move time. [] = estimated automatically

% ---------- Scenario 3: WALL -----------------------------------------------
% A straight wall between the start and finish. The payload has to go around one end.
% To test "stop right next to a wall", put the finish pose close to the wall instead,
% on the same side as the start.
S3.wall_end_A      = [-0.05 0.05];  % m, floor X/Y of one end of the wall
S3.wall_end_B      = [ 0.65 0.05];  % m, floor X/Y of the other end
S3.thickness       = 0.04;          % m, wall thickness
S3.height          = 1.5;           % m, wall height, used for drawing and the clearance check only
S3.start           = [0.8 1.0 -20]; % pose [trolley m, cable m, slew deg]
S3.finish          = [0.4 1.2  30]; % pose
S3.go_around       = 'B';           % which wall end the guess goes around: 'A' or 'B'
S3.guess_clearance = 0.08;          % m, how far past the wall end's keep-out zone the guess goes
S3.tf_guess        = [];            % s, guessed move time. [] = estimated automatically

% ---------- Scenario 4: LIFT WHILE SWINGING --------------------------------
% Slew while raising or lowering the load. The cable length, and so the swing frequency,
% changes during the move, which a fixed-frequency input shaper can't follow.
% The hoist is slow (30 mm/s), so a long lift sets the move time.
S4.start    = [0.6 1.3  0];   % pose [trolley m, cable m, slew deg]
S4.finish   = [0.6 0.7 90];   % pose. Smaller cable = load raised
S4.tf_guess = [];             % s, guessed move time. [] = estimated automatically

% ---------- Scenario 5: TROLLEY + SLEW TOGETHER ----------------------------
% Move the trolley in or out while slewing. Slewing couples the two swing directions,
% which per-axis input shaping ignores.
S5.start       = [0.3  1.0  0];   % pose [trolley m, cable m, slew deg]
S5.finish      = [0.85 1.0 60];   % pose
S5.guess_shape = 'arc';           % starting guess: 'arc' = trolley and slew move together
                                  % (a spiral), 'straight' = payload in a straight line
S5.tf_guess    = [];              % s, guessed move time. [] = estimated automatically

%% ================= Nothing below needs editing ===========================
here = fileparts(mfilename('fullpath'));
addpath(here);
base = struct('name', '', 'start', [], 'finish', [], 'obstacles', struct('kind', {}), ...
              'waypoints', zeros(2, 0), 'guess_shape', 'arc', 'tf_guess', [], ...
              'time_weight', SPEED_WEIGHT, 'preview_only', PREVIEW_ONLY, ...
              'obstacle_stiffness_start', 1e-2, ...
              'out_dir', fullfile(here, 'results'));

switch SCENARIO
    case 1
        n = numel(S1.sweep_angles);
        tf_all = nan(1, n);  swing_all = nan(1, n);
        for i = 1:n
            scn = base;
            ang = S1.sweep_angles(i);
            scn.name     = sprintf('s1_angle_sweep_%gdeg', ang);
            scn.start    = [S1.trolley, S1.cable, S1.start_slew];
            scn.finish   = [S1.trolley, S1.cable, S1.start_slew + ang];
            scn.tf_guess = S1.tf_guess;
            traj = solve_scenario_bvp(scn);
            if ~isempty(traj), tf_all(i) = traj.tf;  swing_all(i) = traj.peak_swing_deg; end
        end
        if ~PREVIEW_ONLY
            T = table(S1.sweep_angles(:), tf_all(:), swing_all(:), ...
                      'VariableNames', {'slew_deg', 'move_time_s', 'peak_swing_deg'});
            fprintf('\nAngle sweep summary:\n');  disp(T);
            writetable(T, fullfile(base.out_dir, 's1_angle_sweep_summary.csv'));
            figure('Name', 'Angle sweep', 'Color', 'w');
            plot(abs(S1.sweep_angles), tf_all, 'o-', 'LineWidth', 1.5);  grid on;
            xlabel('slew angle (deg)');  ylabel('optimal move time (s)');
            title(sprintf('Angle sweep: trolley %.2f m, cable %.2f m', S1.trolley, S1.cable));
        end

    case 2
        scn = base;
        scn.obstacles = pin_obstacle(S2);
        scn.start     = S2.start;
        scn.tf_guess  = S2.tf_guess;
        c = S2.pin_radius + 0.08 + S2.guess_clearance;          % guess radius around the pin center
        switch S2.mode
            case 'pass'
                scn.name   = sprintf('s2_pin_pass_%s', S2.pass_side);
                scn.finish = S2.finish;
                scn.waypoints = pin_pass_waypoints(S2, c);
            case 'loop'
                scn.name   = sprintf('s2_pin_loop_%s', S2.loop_direction);
                scn.finish = S2.start;
                scn.waypoints = pin_loop_waypoints(S2, c);
                scn.obstacle_stiffness_start = 1;   % start == finish: keep the loop from cutting through
            otherwise
                error('S2.mode must be ''pass'' or ''loop''.');
        end
        traj = solve_scenario_bvp(scn);

    case 3
        scn = base;
        scn.name      = sprintf('s3_wall_around_%s', S3.go_around);
        scn.obstacles = struct('kind', 'wall', 'a', S3.wall_end_A, 'b', S3.wall_end_B, ...
                               'thickness', S3.thickness, 'height', S3.height);
        scn.start     = S3.start;
        scn.finish    = S3.finish;
        scn.tf_guess  = S3.tf_guess;
        scn.waypoints = wall_waypoints(S3);
        scn.guess_shape = 'straight';
        traj = solve_scenario_bvp(scn);

    case 4
        scn = base;
        scn.name     = 's4_lift_while_swinging';
        scn.start    = S4.start;
        scn.finish   = S4.finish;
        scn.tf_guess = S4.tf_guess;
        traj = solve_scenario_bvp(scn);

    case 5
        scn = base;
        scn.name        = 's5_trolley_and_slew';
        scn.start       = S5.start;
        scn.finish      = S5.finish;
        scn.guess_shape = S5.guess_shape;
        scn.tf_guess    = S5.tf_guess;
        traj = solve_scenario_bvp(scn);

    otherwise
        error('SCENARIO must be 1 to 5.');
end

if RUN_SIM && ~PREVIEW_ONLY
    run(fullfile(here, 'sim_trajectory.m'));     % plays trajectory.mat (the last solve)
end

%% ---------------- helpers: turn the knobs into a guess path ----------------
function o = pin_obstacle(S2)
o = struct('kind', 'pin', 'xy', S2.pin_xy, 'radius', S2.pin_radius, 'height', S2.pin_height);
end

function P = pose_xy(pose)
P = pose(1) * [cosd(pose(3)); sind(pose(3))];
end

function wp = pin_pass_waypoints(S2, c)
% Three points on a half circle of radius c around the pin, on the chosen side.
pin = S2.pin_xy(:);
radial = pin / norm(pin);                          % from the mast out through the pin
switch S2.pass_side
    case 'outside', side = radial;
    case 'inside',  side = -radial;
    otherwise, error('S2.pass_side must be ''outside'' or ''inside''.');
end
A = pose_xy(S2.start) - pin;  B = pose_xy(S2.finish) - pin;
a_side = atan2(side(2), side(1));
% Directions from the pin toward the start and finish, swept through the chosen side
a0 = atan2(A(2), A(1));  a1 = atan2(B(2), B(1));
d0 = wrap_pi(a_side - a0);  d1 = wrap_pi(a1 - a_side);
ang = [a0 + 0.5*d0, a_side, a_side + 0.5*d1];
wp = pin + c * [cos(ang); sin(ang)];
if norm(pose_xy(S2.start) - pin) < c || norm(pose_xy(S2.finish) - pin) < c
    warning('The start or finish pose is within guess_clearance of the pin; the guess may cut close.');
end
end

function wp = pin_loop_waypoints(S2, c)
% Points evenly spaced around a circle of radius c, starting from the start pose's side.
pin = S2.pin_xy(:);
A = pose_xy(S2.start) - pin;
if norm(A) <= c
    error(['Scenario 2 loop: the start pose is only %.0f mm from the pin center but the loop ' ...
           'radius is %.0f mm. Start farther from the pin or lower guess_clearance.'], 1000*norm(A), 1000*c);
end
switch S2.loop_direction
    case 'ccw', sgn = 1;
    case 'cw',  sgn = -1;
    otherwise, error('S2.loop_direction must be ''ccw'' or ''cw''.');
end
ang = atan2(A(2), A(1)) + sgn * deg2rad(30:60:330);
wp = pin + c * [cos(ang); sin(ang)];
end

function wp = wall_waypoints(S3)
% Two points just past the chosen wall end: one on the start's side, one on the finish's.
a = S3.wall_end_A(:);  b = S3.wall_end_B(:);
switch S3.go_around
    case 'A', E = a;  d = (a - b) / norm(a - b);
    case 'B', E = b;  d = (b - a) / norm(b - a);
    otherwise, error('S3.go_around must be ''A'' or ''B''.');
end
n = [-d(2); d(1)];                                  % normal to the wall
c = S3.thickness/2 + 0.08 + S3.guess_clearance;
sA = sign(n.' * (pose_xy(S3.start)  - E));
sB = sign(n.' * (pose_xy(S3.finish) - E));
if sA == sB
    wp = zeros(2, 0);                               % same side: no need to go around
    fprintf('Scenario 3: start and finish are on the same side of the wall, so the guess goes straight.\n');
else
    wp = [E + c*d + c*sA*n, E + c*d + c*sB*n];
end
end

function a = wrap_pi(a)
a = mod(a + pi, 2*pi) - pi;
end
