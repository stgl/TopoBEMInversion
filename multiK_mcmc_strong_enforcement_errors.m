filename = 'multiK_mcmc_strong_enforcement_1Msamples';
startIndex = 500000;
numberOfRandomSamples = 3000;

load(filename);
samples = x_keep(startIndex:end,:);
subsamples = datasample(samples,numberOfRandomSamples,'Replace',false);

% Setting the correct directories
p=pathdef; path(p)

% Concavity:
theta = 0.4;

% weighting factors for bay and Loma Prieta points:
w_bay_constr = 10^(-6.5);
w_lp_constr = 10^(-6.8);

sig_elev=16; % in m

fission_track_uplift=0.0008;   %m/yr, @ Loma Prieta from Roland's rise
                                % and fall paper

multiple_K_flag = 1; % 0 = single K for every point, 1 = use lithology.

minLogK = -8;
maxLogK = -2;

Niter = 1E5;
xstep = 1E-2;

% NO CHANGES BELOW HERE.
[geo_map, e_outlets, Ginv_elev, channel_indexes, e_chan, n_K, ...
    ind_chan_misfit, litho_chan, outlet_indexes] = ...
    inputs_for_objective_function(theta);

G_chan = load('input/G_tect_shear_push.txt');
G_bay = load('input/G_bay_shear_push.txt');
G_lp = load('input/G_lomaprieta_shear_push.txt');

% Take account of the 2 BC due to 1 imposed on each side of the box
G_chan=G_chan./2; G_bay=G_bay./2; G_lp=G_lp./2;

bay_constr=zeros(length(G_bay(:,1)));
lp_constr=ones(length(G_lp(:,1)),1).*fission_track_uplift;

fun = @(K)topo_linear_lsq_model_cov(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

v = zeros(length(subsamples), 5);
for i = 1:length(subsamples)
  [thisv, thiscov_v] = fun(subsamples{i});
  v(i,1:2) = thisv;
  v(i,3) = thiscov_v(1,1);
  v(i,4) = thiscov_v(2,2);
  v(i,5) = thiscov_v(1,2);
end

save multiK_mcmc_strong_enforcement_v_with_cov v
