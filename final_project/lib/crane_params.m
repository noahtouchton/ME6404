function p = crane_params(cfg)
% Parameter vector for the generated files, in the exact order of params_sym in EOM.m.
% Building it from the cfg struct (instead of an anonymous function) means the saved
% settings and the numbers the solver used can never disagree.
p = [cfg.g; cfg.r1; cfg.r2; cfg.r3; cfg.xmax; cfg.lmax; cfg.dxmax; cfg.dlmax; cfg.dthmax; ...
     cfg.ddxmax; cfg.ddlmax; cfg.ddthmax; cfg.R_p; 0; cfg.r_p; cfg.c_clr; ...
     cfg.W_obs; cfg.k_x; cfg.W_t; cfg.k_v];
end
