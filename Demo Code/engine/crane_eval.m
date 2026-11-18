function Z = crane_eval(S, t)
% Evaluate a saved solution S at real times t (22 x numel(t)). Solutions live on tau = t/t_f.
Z = deval(S.sol, min(max(t / S.t_f, 0), 1));
end
