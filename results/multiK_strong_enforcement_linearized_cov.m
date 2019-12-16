scenario = 'strong';
n = 5000;

% Setting the correct directories
p=pathdef; path(p)

load(strcat('results/multiK_', scenario, '_enforcement'));

% Calcualte Jacobian:

J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

covK = J * sig_elev.^2;

% Draw random samples:
samples = mvnrnd(K,covK,n);

fun = @(K)topo_linear_lsq_model_cov(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);
v = zeros(length(samples), 5);
for i = 1:length(samples(:,1))
  [thisv, thiscov_v] = fun(ssamples(i,:));
  v(i,1:2) = thisv;
  v(i,3) = thiscov_v(1,1);
  v(i,4) = thiscov_v(2,2);
  v(i,5) = thiscov_v(1,2);
end

save strcat('results/multiK_', scenario, '_enforcement_v_with_cov_linearized') v
