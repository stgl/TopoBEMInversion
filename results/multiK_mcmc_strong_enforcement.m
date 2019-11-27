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

% set Ko equal to fminsearch solution for best fit (this will be adjusted to make sure that all values are within prior bounda

Ko = [-4.8494   -4.5998   -4.7279    0.7075   -4.5056   -5.0411   -5.1443 -4.5885   -4.7222   -4.4796   -4.5087   -2.8392   -4.8032   -4.6731 -4.4013   -4.6969   -5.0556   -4.5692   -5.0758   -1.4962   -4.4038 -4.5888   -4.1285   -4.1497   -2.1194   -4.6673   -4.7326   -1.5071 -1.9293];


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

n_constraints = 2;
n_Gtect_rows = length(G_chan(:,1));
n_Gtect_cols = length(G_chan(1,:));
n_outlets = length(Ginv_elev(1,:))-n_Gtect_rows-n_constraints;
n_rows = n_Gtect_rows + n_outlets + n_constraints;
n_cols = n_Gtect_cols + n_outlets;

% Set D for MCMC:

D.d = [e_chan;e_outlets;bay_constr;lp_constr];

invSig = sparse(ind_chan_misfit,ind_chan_misfit, ...
    ones(1,length(ind_chan_misfit)),n_rows,n_rows) .* sig_elev.^2;
invSig(n_Gtect_rows+1:n_Gtect_rows+n_outlets, ...
    n_Gtect_rows+1:n_Gtect_rows+n_outlets) = ...
    speye(n_outlets, n_outlets) * sig_elev.^2;
invSig(n_Gtect_rows+n_outlets+1,n_Gtect_rows+n_outlets+1) =  w_bay_constr.^2;
invSig(n_rows, n_rows) = w_lp_constr.^2;

D.invSig = invSig;

% Set X for MCMC:

Ko = (Ko > maxLogK).*(maxLogK - xstep) + (Ko < minLogK).*(minLogK + xstep) + ( (Ko <= maxLogK) & (Ko >= minLogK) ).*Ko;

Ko = Ko';

X.x0     = Ko;
X.xstep  = [xstep*ones(length(Ko),1)];
X.xbnds  = [minLogK*ones(length(Ko),1) maxLogK*ones(length(Ko),1)];
X.xprior = [];
X.C = Inf; % TODO: This needs to be checked by Andreas!!!

tic

fun = @(K)topo_linear_lsq_wrapper(K, e_chan, e_outlets, ...
    sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
    w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

[x_keep, logL_keep, logQ_keep, accrate] = mcmc(fun,D,X,Niter,1);

t=toc;

disp(['Running MCMC routine takes ',num2str(t),' secs'])

save multiK_mcmc_strong_enforcement
quit
