function [v, cov_v] = topo_linear_lsq_model_cov(K, e_chan, ...
e_outlets, sig_elev, ind_chan_misfit, G_chan, Ginv_elev, ...
bay_constr, w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map, varargin)

if(length(varargin) == 0)
  n = 1.0;
elseif(length(varargin) == 1)
  n = varargin{1};
else
  error('Improper number of arguments to topo_linear_lsq_soln.m');
end

d = [e_chan;e_outlets;bay_constr.^(1./n);lp_constr.^(1./n)];

logKcon = K;

% Expand G_tect matrix to account for outlets:

n_constraints = 2;
n_Gtect_rows = length(G_chan(:,1));
n_Gtect_cols = length(G_chan(1,:));
n_outlets = length(Ginv_elev(1,:))-n_Gtect_rows-n_constraints;
n_rows = n_Gtect_rows + n_outlets + n_constraints;
n_cols = n_Gtect_cols + n_outlets;

Gt = sparse(n_rows,n_cols);
Gt(1:n_Gtect_rows,1:n_Gtect_cols) = G_chan.^(1./n);
Gt(n_Gtect_rows+1:n_Gtect_rows+n_outlets, ...
    n_Gtect_cols+1:n_Gtect_cols+n_outlets) = speye(n_outlets,n_outlets);
Gt(n_Gtect_rows+n_outlets+1:n_Gtect_rows+n_outlets+n_constraints, ...
    1:n_Gtect_cols) = [G_bay.^(1./n);G_lp.^(1./n)];



% Greate K matrix to account for steepness:

K = sparse(n_rows, n_rows);

if(length(logKcon) == 1)
    K(1:n_Gtect_rows,1:n_Gtect_rows) = ...
        speye(n_Gtect_rows, n_Gtect_rows) * (1 ./ 10.^logKcon(1));
else
    for(i=1:length(geo_map(:,1)))
        n_el = geo_map(i,3) - geo_map(i,2) + 1;
        K(geo_map(i,2):geo_map(i,3),geo_map(i,2):geo_map(i,3)) = ...
            speye(n_el,n_el) * (1 ./ 10.^logKcon(i));
    end
end
K(n_Gtect_rows+1:n_Gtect_rows+n_outlets, ...
    n_Gtect_rows+1:n_Gtect_rows+n_outlets) = speye(n_outlets,n_outlets);
K(n_Gtect_rows+n_outlets+1:n_rows, ...
    n_Gtect_rows+n_outlets+1:n_rows) = speye(n_constraints);

K = K^(1./n);

% Create Weighting Matrix:

W = sparse(ind_chan_misfit,ind_chan_misfit, ...
    ones(1,length(ind_chan_misfit)),n_rows,n_rows) / sig_elev.^2;
W(n_Gtect_rows+1:n_Gtect_rows+n_outlets, ...
    n_Gtect_rows+1:n_Gtect_rows+n_outlets) = ...
    speye(n_outlets, n_outlets) / sig_elev.^2;
W(n_Gtect_rows+n_outlets+1,n_Gtect_rows+n_outlets+1) = 1 / w_bay_constr.^(2/n);
W(n_rows, n_rows) = 1 / w_lp_constr.^(2/n);

% Create Data Covariance Matrix:

cov_d = sparse(ind_chan_misfit,ind_chan_misfit, ...
    ones(1,length(ind_chan_misfit)),n_rows,n_rows) * sig_elev.^2;
cov_d(n_Gtect_rows+1:n_Gtect_rows+n_outlets, ...
    n_Gtect_rows+1:n_Gtect_rows+n_outlets) = ...
    speye(n_outlets, n_outlets) * sig_elev.^2;
cov_d(n_Gtect_rows+n_outlets+1,n_Gtect_rows+n_outlets+1) = w_bay_constr.^(2/n);
cov_d(n_rows, n_rows) = w_lp_constr.^(2/n);

G = Ginv_elev*K*Gt;

lastwarn('')

m = (G'*W*G)\G'*W*d;

cov_m = (G'*G)\G'*cov_d*((G'*G)\G')';

v = m(1:2).^n;
cov_v = cov_m(1:2,1:2).^n;
