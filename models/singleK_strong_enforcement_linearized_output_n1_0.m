scenario = 'strong';
n = 1.0;

% Setting the correct directories
p=pathdef; path(p)

load(strcat('results/singleK_', scenario, '_enforcement_n', num2str(n)));

% Calcualte Jacobian:

J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, n);

covK = J * sig_elev.^2;

for(i=1:length(K))
  fprintf('K %i (in order): %6.2f ± %6.4f\n', i, K(i), sqrt(covK(i,i)));
end
