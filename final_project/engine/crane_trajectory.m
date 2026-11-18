function T = crane_trajectory(S, t)
% States, optimal inputs u*, Hamiltonian and geometry of saved solution S at real times t.
p = crane_params(S.cfg);
Z = crane_eval(S, t);
N = numel(t);
T.t = t;
T.X = Z(1:11, :);
T.U = zeros(3, N);
T.H = zeros(1, N);
for k = 1:N
    T.U(:, k) = optimal_inputs(Z(:, k), p);
    T.H(k)    = hamiltonian_value(Z(:, k), p);
end
T.G = crane_geometry(T.X, S.cfg);
end
