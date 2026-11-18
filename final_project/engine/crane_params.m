function p = crane_params(cfg)
% Parameter vector for the generated functions. The order comes from
% generated/param_names.m (written by theory/EOM.m), so it can't drift out of sync.
n = param_names();
p = zeros(numel(n), 1);
for k = 1:numel(n)
    p(k) = cfg.(n{k});
end
end
