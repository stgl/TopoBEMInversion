scenario = 'strong';
n = 5000;
max_range = 5;

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
max_samples = repmat(K,n,1) + max_range / 2;
min_samples = repmat(K,n,1) - max_range / 2;

i = find(samples > max_samples);
samples(i) = max_samples(i);
i = find(samples < min_samples);
samples(i) = min_samples(i);

fun = @(K)topo_linear_lsq_model_cov(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);
v = zeros(length(samples), 5);
parfor i = 1:length(samples(:,1))
  [thisv, thiscov_v] = fun(samples(i,:));
  v(i,1:2) = thisv;
  v(i,3) = thiscov_v(1,1);
  v(i,4) = thiscov_v(2,2);
  v(i,5) = thiscov_v(1,2);
  if(mod(i,100))
    fprintf('%i...', i);
  end
end

filename = strcat('results/multiK_', scenario, '_enforcement_v_with_cov_linearized');
save(filename, 'v');
