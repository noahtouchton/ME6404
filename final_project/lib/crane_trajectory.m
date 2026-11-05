function T = crane_trajectory(S, t)
% States, optimal inputs, Hamiltonian and geometry of saved solution S at times t.
params = crane_params(S.cfg);
Z = crane_eval(S, t);
N = numel(t);

T.t = t;
T.X = Z(1:11, :);
T.U = zeros(3, N);
T.H = zeros(1, N);
for k = 1:N
    T.U(:, k) = optimal_inputs(Z(:, k), params);
    T.H(k)    = hamiltonian_value(Z(:, k), params);
end
T.G = crane_geometry(T.X, S.cfg);
end
