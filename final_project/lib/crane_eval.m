function Z = crane_eval(S, t)
% Evaluate a saved solution S at real times t (22 x numel(t)).
% Free-tf solutions live on tau = t/t_f in [0, 1]; fixed-tf ones on real time.
if S.time_scaled
    Z = deval(S.sol, t / S.t_f);
else
    Z = deval(S.sol, t);
end
end
