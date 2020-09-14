n = 0.666667;
scenario = 'weak';

% Setting the correct directories
p=pathdef; path(p)

% Concavity:
theta = 0.4;

% weighting factors for bay and Loma Prieta points:
w_bay_constr = 10^(-6.3);
w_lp_constr = 10^(-6.1);

sig_elev=16; % in m

fission_track_uplift=0.0008;   %m/yr, @ Loma Prieta from Roland's rise
                                % and fall paper

multiple_K_flag = 0; % 0 = single K for every point, 1 = use lithology.

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

if(multiple_K_flag == 1)
    Ko = -4 .* ones(1,length(geo_map(:,1)));
else
    Ko = -4;
end

opts = optimset('MaxIter', 1E6, 'MaxFunEvals', 1E6);

tic

fun = @(K)topo_linear_lsq_misfit(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, n);
[K, fval] = fminsearch(fun,Ko, opts);

t=toc;

[x_channels, y_channels]=getCoordinates;

x_channels=x_channels(channel_indexes);
y_channels=y_channels(channel_indexes);

[log10_wrss_elev_outlets,elev_pred_outletsfree, outlets_pred, vshear, ...
    vconverge, c_bay, c_lp] = topo_linear_lsq_stats(K, e_chan, ...
    e_outlets, sig_elev, ind_chan_misfit, G_chan, Ginv_elev, ...
    bay_constr, w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, n);

outlets_resid = e_outlets - outlets_pred;

residuals_elev_outletsfree=e_chan-elev_pred_outletsfree;

nr_data=length(ind_chan_misfit);
nr_model_par = length(K)+length(e_outlets)+2;

chi_sq_outletsfree=(1/((nr_data-nr_model_par)*sig_elev^2))*...
    (sum((e_chan-elev_pred_outletsfree).^2));


fprintf('Chi-2 misfit: %6.2f\n', chi_sq_outletsfree);
fprintf('Shear displacement rate (mm/yr): %6.2f\n', vshear*1000);
fprintf('Convergent displacement rate (mm/yr): %6.2f\n', vconverge*1000);
fprintf('Uplift rate at Bay constraint (mm/yr): %6.2f\n', c_bay*1000);
fprintf('Uplift rate at LP (mm/yr): %6.2f\n', c_lp*1000);
fprintf('Log10 K value (log m/yr): %6.2f\n', K);

filename = strcat('singleK_', scenario, '_enforcement_n',num2str(n));

save(filename)

quit
