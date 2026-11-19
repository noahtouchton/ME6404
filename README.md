# ME6404

## Final project: optimal swing-around maneuver for the tower crane

The payload has to get around a pin that the trolley can't reach, so it has to swing around
it, then return to rest at its starting state. Two methods make that trajectory:

- **Optimal control** (Pontryagin's minimum principle) solves for the whole loop at once: the
  states, their costates and the input form a two-point boundary-value problem, solved with
  `bvp4c` within the hardware limits on trolley, hoist and slew speed and acceleration.
- **Input shaping** builds the loop from three point-to-point moves and plays it forward. Each
  move is shaped so the swing it starts is cancelled again by the end.

### Layout

```
final_project/
├── crane.m              shared helpers, as one class of static functions: the task, model
│                        parameters, shapers and shaped moves, simulation, trajectories, checks
├── solve_oc.m           optimal control: fixed 8 s loop, then minimum time
│                        -> results/oc_fixed.mat, oc_min_time.mat
├── solve_is.m           input shaping: classic, cancelling and robust loops -> results/is_*.mat
├── compare_methods.m    replays every saved loop, as designed and with the wrong cable length
│                        -> results/comparison.csv, comparison.mat, results/figures/*.png
├── run_sim.m            one saved loop: checks, plots, 5 s coast, 3D animation (no solving)
├── lab5_trajectories.m  input-shaping test files for the Lab 5 tower crane -> lab5/
├── model/
│   ├── EOM.m            symbolic dynamics, costates and optimal input -> writes model/generated/
│   └── generated/       written by EOM.m (don't edit by hand)
├── lab5/                Lab 5 upload files (.xlsx), summary.csv, predicted.png
└── results/             saved loops, comparison table and figures (old/: earlier versions)
```

A saved loop (`results/*.mat`) holds the task (`cfg`), its kind (`optimal` or `commands`), the
final time, a label, and either the `bvp4c` solution or the shaped command: a list of
acceleration ramps, one row `[axis, start, duration, acceleration]` each.

### Input-shaped loop

`solve_is.m` builds the loop from the same three moves every time. Each move is the fastest
point-to-point profile (velocity steps made as ramps at the acceleration limit) convolved with a
shaper:

1. **stage**: slew to -th_s while the trolley runs out to `xmax`;
2. **sweep**: slew from -th_s to +th_s at 75-100 % of the slew limit (scanned) with the trolley at
   `xmax`. Slewing at full speed throws the payload outward (centrifugal, about 31 mm) past
   `xmax`, and clearing the pin on the outside takes 25 mm;
3. **return**: trolley back in and slew back to 0, passing inside the pin.

Shaping each axis on its own ignores how slewing couples them (the trolley's ground position is
x·e^(iθ)), so the designs differ in how they deal with the swing that coupling leaves:

- **classic**: standard ZV, ZVD or EI shapers; shaper type, sweep speed and sweep length come from
  a scan on a fast linear model. Some swing is left.
- **cancel**: the stage and return shapers are solved on the full nonlinear model so the payload is
  at rest after the stage and at the end (two impulses per axis, zero swing at the design cable
  length).
- **robust**: like cancel, with three-impulse shapers fitted over cable lengths l0 ± 5 % (ZVD-like).

No Optimization Toolbox is needed: the shapers are solved with a small Levenberg-Marquardt routine
inside `solve_is.m`.

### Results

From `compare_methods`. Every command is replayed open loop on the full nonlinear model, and all
were designed for l0 = 1.00 m. Margins are distances beyond the required clearance. "Band" is the
range of real cable lengths, on a 2 % grid, over which the replay stays admissible and leaves
< 5 mm of swing.

| Method | t_f (s) | Past xmax (mm) | Pin margin (mm) | Cable margin (mm) | Swing left (mm) | Band, l (m) |
|---|---:|---:|---:|---:|---:|:---:|
| OC fixed 8 s | 8.00 | 39.5 | 14.5 | 9.4 | 0.0 | 1.00-1.00 |
| OC min-time | 6.85 | 32.5 | 2.4 | 4.4 | 0.0 | 1.00-1.00 |
| IS classic (EI) | 19.03 | 44.2 | 10.3 | 21.6 | 7.8 | none |
| IS cancel (ZV) | 13.28 | 41.2 | 10.2 | 21.4 | 0.0 | 0.98-1.00 |
| IS cancel robust (ZVD) | 18.83 | 35.8 | 10.4 | 21.6 | 0.7 | 0.88-1.06 |

With the wrong cable length (same commands):

| Method | Swing left (mm) at l = 0.90 / 0.96 / 1.00 / 1.04 / 1.10 m | Pin margin (mm) at the same l | Worst pin margin, 0.8-1.2 m |
|---|---|---|---:|
| OC fixed 8 s | 39.5 / 17.7 / 0.0 / 19.8 / 52.7 | -5.6 / 7.3 / 14.5 / 19.2 / 19.7 | -28.8 |
| OC min-time | 105.6 / 44.5 / 0.0 / 45.3 / 111.5 | -41.3 / -16.3 / 2.4 / -9.3 / -24.7 | -48.7 |
| IS classic (EI) | 18.6 / 11.6 / 7.8 / 15.5 / 32.5 | 2.3 / 6.0 / 10.3 / 15.9 / 25.5 | 1.6 |
| IS cancel (ZV) | 24.6 / 8.4 / 0.0 / 12.0 / 38.1 | 8.0 / 9.2 / 10.2 / 11.7 / 14.6 | 3.9 |
| IS cancel robust (ZVD) | 2.7 / 1.9 / 0.9 / 1.9 / 7.0 | 5.9 / 8.6 / 10.4 / 12.4 / 15.9 | 3.0 |

What it shows:

- **Time:** optimal control is fastest, 6.85 s at minimum time against 13.28 s for the fastest
  shaped loop here.
- **How:** the minimum-time optimal loop keeps the trolley within 0.85 m and pumps up a swing of up
  to 159 mm; the shaped loops drive the trolley to its 0.925 m limit and let the slew throw the
  payload out.
- **Robustness:** off the design cable length, the robust shaped loop keeps the leftover swing at
  or below 5 mm from 0.88 to 1.06 m. Every other method leaves 8-45 mm at 0.96 or 1.04 m and up to
  112 mm at 0.90 or 1.10 m.
- **Pin margin:** the minimum-time optimal loop rides close to the clearance (2.4 mm) and clips it
  for any cable-length error of 4 % or more. The shaped loops keep a positive margin from 0.8 to 1.2 m, because their
  reach comes from slewing and a longer cable swings farther out.
- **Classic shaping is not enough:** one standard shaper per axis ignores slewing's coupling and
  leaves swing. Shapers solved for the coupled crane fix it.

Figures: `results/figures/paths.png`, `swing.png`, `robustness.png`, `actuators.png`.

### Lab 5 test trajectories

`lab5_trajectories.m` writes continuous-velocity upload files for the Lab 5 tower crane
(`lab5/*.xlsx`: first sheet, headers `time_s, trolley_pct, slew_pct, hoist_pct` in row 1, one row
every 0.05 s, velocity in % of the Table IV maximum, ending at zero) that test whether input
shaping removes the swing on the real crane. Each row is the average velocity over its step, so the
positions add up exactly. Each test starts from the Lab 5 start pose (trolley 550 mm, slew 195°,
cable 1350 mm) and makes one move.

The lab crane is a double pendulum (0.357 kg hook on 1.35 m of cable, 0.050 kg payload 0.381 m
below it), with swing periods of about 2.38 s and 1.13 s, so the shapers target both modes
(`ZV12`: ZV for each mode, convolved). All times sit on the 0.05 s grid. Predictions come from a
linear double-pendulum model replaying each file:

| File | Shaper | Speed | Move | Predicted swing left (peak, mm) | Predicted translation penalty (s) |
|---|---|---:|---|---:|---:|
| trolley_unshaped | none | 100 % | 300 mm out | 38.8 | 7.14 |
| trolley_ZV_mode1 | ZV, mode 1 | 100 % | 300 mm out | 1.3 | 0.18 |
| trolley_ZV_both | ZV12 | 100 % | 300 mm out | 0.4 | 0.06 |
| slew_unshaped | none | 100 % | +40° | 223.9 | 47.05 |
| slew_ZV_both | ZV12 | 100 % | +40° | 32.2 | 6.16 |
| slew_ZV_both_half | ZV12 | 50 % | +40° | 7.2 | 1.67 |

Things to know before running them:

- Run `trolley_unshaped` first, read the two swing periods off the recorded data, put them in
  `T_meas` at the top of the script, and regenerate. The shapers are only as good as those periods.
- The slew tests turn away from the finish (195° to 235°) to stay clear of the T obstacle; check
  that side is clear. `slew_unshaped` swings the payload about ±22 cm on the model.
- Shaping the slew can't remove the outward throw from slewing itself (it grows with speed
  squared), which is why the full-speed shaped slew still leaves about 32 mm and half speed about
  7 mm.
- The predictions assume the crane holds each row's velocity for its 0.05 s step; check the
  interface's "File formats & limits" panel.

### Usage

Run from `final_project/` in MATLAB (Symbolic Math Toolbox needed only for `model/EOM.m`):

1. `model/EOM.m` once, or after changing the model or cost: writes `model/generated/`.
2. `solve_oc` (about 5 minutes): fixed-time loop, then minimum time.
3. `solve_is` (about 7 minutes): classic, cancelling and robust shaped loops.
4. `compare_methods`: the table and figures above.
5. `run_sim`: look at one loop (set `solution` at the top).
6. `lab5_trajectories`: the Lab 5 test files.
