function res = bc_loop(ya, yb, X0)
% 22 boundary conditions: start at X0, eta = 0; end at X0; p_eta = 0.
% eta(tf) = 2*pi*kdir is NOT imposed: eta is a function of the path's endpoints
% and topology, so that BC is degenerate and leaves p_eta undetermined (singular
% Jacobian). The loop count is set by the guess + the pin wall; the solve
% wrappers check eta(tf) afterwards.
    res = [ya(1:10) - X0;
           ya(11);
           yb(1:10) - X0;
           yb(22)];
end
