function J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, varargin)

if(length(varargin) == 0)
  n = 1.0;
elseif(length(varargin) == 1)
  n = varargin{1};
else
  error('Improper number of arguments to topo_linear_lsq_soln.m');
end

log_dK = 1E-10;
inv_J = zeros(length(K),length(e_chan));

% Calculate base d:

d = [e_chan;e_outlets;bay_constr.^(1./n);lp_constr.^(1./n)];

[~, dpred_center, W] = topo_linear_lsq_soln(d, K, sig_elev, ind_chan_misfit, ...
    G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map, n);

z_center = dpred_center(1:length(e_chan));

for(i=1:length(K))
  thisK = K;
  thisK(i) = thisK(i) + log_dK;
  [~, this_dpred, W] = topo_linear_lsq_soln(d, thisK, sig_elev, ind_chan_misfit, ...
      G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map, n);
  z_this = this_dpred(1:length(e_chan));
  dzdK = (z_this - z_center) / log_dK;
  inv_J(i,:) = dzdK';
end

J = inv(inv_J*inv_J');
