function J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map)

log_dK = 1E-6;
J = zeros(length(K),length(e_chan));

% Calculate base d:

d = [e_chan;e_outlets;bay_constr;lp_constr];

[~, dpred_center, W] = topo_linear_lsq_soln(d, K, sig_elev, ind_chan_misfit, ...
    G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map);

z_center = dpred_center(1:length(e_chan));

for(i=1:length(K))
  thisK = K;
  thisK(i) = thisK(i) + log_dK;
  [~, this_dpred, W] = topo_linear_lsq_soln(d, K, sig_elev, ind_chan_misfit, ...
      G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map);
  z_this = this_dpred(1:length(e_chan));
  dKdz = ones(length(z_this),1)* log_dK ./ (z_this - z_center);
  J(i,:) = dKdz';
end
