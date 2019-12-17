scenario = 'strong';

% Setting the correct directories
p=pathdef; path(p)

load(strcat('results/multiK_', scenario, '_enforcement'));

% Calcualte Jacobian:

J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

covK = J * sig_elev.^2;

for(i=1:length(K))
  fprintf('K %i (in order): %6.2f ± %6.2f', i, K(i), covK(i,i));
end
