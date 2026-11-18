# ME6404

## Final project: optimal swing-around maneuver for the tower crane

The payload has to get around a pin that the trolley can't reach around, so it
**swings** behind it, then settles back at its start with no residual
vibration. The trajectory comes straight from **Pontryagin's minimum principle**:
a two-point boundary value problem solved with `bvp4c` under the fixed hardware
limits. It's exported as the crane GUI's continuous-velocity Excel upload
(`time_s, trolley_pct, slew_pct, hoist_pct` at 0.05 s).

### Layout

```
final_project/
├── lab_tune.m          ← the lab script: maneuver + lab corrections -> solve -> check -> .xlsx
├── run_sim.m           ← replay/animate the last solution (no solving)
├── theory/             ← the math, readable
│   ├── EOM.m             dynamics -> cost -> H -> u* (∂H/∂u = 0) -> costates (λ' = -∂H/∂x)
│   └── pmp_bvp.m         boundary conditions + transversality (H(t_f) = 0) handed to bvp4c
├── engine/             ← behind the scenes: things you shouldn't need to touch
│   ├── hidden_settings.m hardware limits, cost weights, wall stiffness, solver tolerances
│   ├── solve_maneuver.m  warm/cold start, continuation, wall hardening, minimum time
│   └── ...               problem setup, initial guess, geometry, checks, preview plot
├── export/             ← the .xlsx upload + an open-loop replay of exactly that file
├── generated/          ← written by theory/EOM.m (don't edit)
├── results/            solution.mat: the last solve, used as the warm start (git-ignored)
└── uploads/            the .xlsx files (+ the .mat that made each one) (git-ignored)
```

### In the lab

1. Open `final_project/` in MATLAB and run `lab_tune` (cold start ~10 s, warm re-solve ~2 s).
2. Upload the printed `uploads/trajectory_*.xlsx`.
3. Watch the run. Change the knob in section 2 of `lab_tune.m` that matches what
   you saw (the symptom → knob table is in the file), then run it again. If only
   `axis_gain` changed, it re-exports without re-solving.

After changing the model or the cost in `theory/EOM.m`, run it once to regenerate
`generated/`. The parameter order is shared automatically through
`generated/param_names.m`.
