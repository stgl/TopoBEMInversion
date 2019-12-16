function J = topo_jacobian_K(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map)

dE = 1E-6;
J = zeros(length(K),length(e_chan));

base_fun = @(Ki, e_chan)topo_linear_lsq_misfit(Ki, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

Ko = K;

for(i=1:length(e_chan))

  fprintf('%i/%i',i,length(e_chan));
  this_e_chan = e_chan;
  this_e_chan(i) = this_e_chan(i) + dE;
  this_fun = @(Kt)base_fun(Kt, this_e_chan);
  [this_K, fval] = fminsearch(this_fun, Ko, opts);

  dKdz = (this_K - Ko) / dE;
  J(i,:) = dKdz;

end
