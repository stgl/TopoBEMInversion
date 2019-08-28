% Concavity:
theta = 0.4;

save_flag=1;

% weighting factors for bay and Loma Prieta points:
log_w_bay = -7.7:0.1:-5.7;
log_w_lp = -7:0.1:-5;

w_bay_constr = 10^(-6.5);
w_lp_constr = 10^(-6.5);

sig_elev=16; % in m

fission_track_uplift=0.0008;   %m/yr, @ Loma Prieta from Roland's rise 
                                % and fall paper

multiple_K_flag = 0; % 0 = single K for every point, 1 = use lithology.
                                
% NO CHANGES BELOW HERE.

[Log_w_bay, Log_w_lp] = meshgrid(log_w_bay, log_w_lp);


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

[x_channels, y_channels]=getCoordinates;

x_channels=x_channels(channel_indexes); 
y_channels=y_channels(channel_indexes);

chi_sq_outletsfree =NaN(size(Log_w_bay));
v_s = NaN(size(Log_w_bay));
v_n = NaN(size(Log_w_bay));
k = NaN(size(Log_w_bay));
bay_resid = NaN(size(Log_w_bay));
lp_resid = NaN(size(Log_w_bay));
t_tot=0;
total_iter=length(log_w_bay)*length(log_w_lp);
nr_iter=0;

for i=1:length(Log_w_bay(:,1))
    for j=1:length(Log_w_bay(1,:))
        tic
        w_bay_constr = 10.^Log_w_bay(i,j);
        w_lp_constr = 10.^Log_w_lp(i,j);
        fun = @(K)topo_linear_lsq_misfit(K, e_chan, e_outlets, ...
            sig_elev, ind_chan_misfit, G_chan, Ginv_elev, bay_constr, ...
            w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);
        try
            [K, fval] = fminsearch(fun,Ko);

            [log10_wrss_elev_outlets,elev_pred_outletsfree, outlets_pred, vshear, ...
                vconverge, c_bay, c_lp] = topo_linear_lsq_stats(K, e_chan, ...
                e_outlets, sig_elev, ind_chan_misfit, G_chan, Ginv_elev, ...
                bay_constr, w_bay_constr, G_bay, lp_constr, w_lp_constr, G_lp, geo_map);

            outlets_resid = e_outlets - outlets_pred;

            residuals_elev_outletsfree=e_chan-elev_pred_outletsfree;

            nr_data=length(ind_chan_misfit);
            nr_model_par = length(K)+length(e_outlets)+2;

            chi_sq_outletsfree(i,j)=(1/((nr_data-nr_model_par)*sig_elev^2))*...
                (sum((e_chan-elev_pred_outletsfree).^2));
            v_s(i,j) = vshear;
            v_n(i,j) = vconverge;
            k(i,j) = K(1);
        catch
            chi_sq_outletsfree(i,j) = nan;
            v_s(i,j) = nan;
            v_n(i,j) = nan;
            k(i,j) = nan;
        end
 
        fprintf('Done with Bay: %6.2f, LP: %6.2f\n', Log_w_bay(i,j), Log_w_lp(i,j))
        
        bay_resid(i,j)=bay_constr-c_bay;
        lp_resid(i,j)=lp_constr-c_lp;
        
        nr_iter=nr_iter+1;
        t=toc;
        t_tot=t_tot+t;
        
        fprintf('Done with Bay: %6.2f, LP: %6.2f\n',...
            Log_w_bay(i,j), Log_w_lp(i,j))
        fprintf('Elapsed time: %6.3f minutes, iter: %d out of %d\n',...
            t_tot/60, nr_iter, total_iter)
    end
end


save(['WeightingResults',num2str(save_flag),'.mat'],'log_w_bay','log_w_lp',...
    'chi_sq_outletsfree','v_s','v_n','k','bay_resid','lp_resid','t_tot');
fprintf('Total time: %6.2f minutes\n', t_tot/60)

