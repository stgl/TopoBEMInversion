function [log_wrss_elev, elev_pred, outlet_pred, vshear, vconverge, ...
    bay_constraint, lp_constraint] = topo_linear_lsq_stats(K, e_chan, ...
    e_outlets, sig_elev, ind_chan_misfit, G_chan, Ginv_elev, ...
    bay_constr, w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map)

d = [e_chan;e_outlets;bay_constr;lp_constr];

[m, dpred, W] = topo_linear_lsq_soln(d, K, sig_elev, ind_chan_misfit, ...
    G_chan, Ginv_elev, w_bay_constr, G_bay, w_lp_constr, G_lp, geo_map);

n_constraints = 2;
n_Gtect_rows = length(G_chan(:,1));
n_outlets = length(Ginv_elev(1,:))-n_Gtect_rows-n_constraints;
n_rows = n_Gtect_rows + n_outlets + n_constraints;

Wmod = W;
Wmod(n_Gtect_rows+n_outlets+1:n_rows,n_Gtect_rows+n_outlets+1:n_rows) = 0.0;

log_wrss_elev=log10(sum(Wmod*(d-dpred).^2));
elev_pred = dpred(1:n_Gtect_rows);
outlet_pred = dpred(n_Gtect_rows+1:n_Gtect_rows+n_outlets);
vshear = m(1);
vconverge = m(2);
bay_constraint = dpred(end-1);
lp_constraint = dpred(end);
