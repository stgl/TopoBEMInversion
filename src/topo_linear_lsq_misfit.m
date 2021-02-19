function log_wrss_comb=topo_linear_lsq_misfit(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, varargin)

if(length(varargin) == 0)
  n = 1.0;
elseif(length(varargin) == 1)
  n = varargin{1};
else
  error('Improper number of arguments to topo_linear_lsq_soln.m');
end

d = [e_chan;e_outlets;bay_constr;lp_constr];

[~, dpred, W] = topo_linear_lsq_soln(d, K, sig_elev, ind_chan_misfit, ...
    G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map, n);

log_wrss_comb=log10(sum(W*(d-dpred).^2));
