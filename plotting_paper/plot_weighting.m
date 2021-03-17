addpath(genpath('/Users/felipearon/Dropbox/Documents/Codes/colorpalettes'))
addpath(genpath('/Users/felipearon/Dropbox/Documents/Codes/utils'))
addpath(genpath('/Users/felipearon/Dropbox/Documents/GitHub/TopoBEMInversion'))

% filename_weight='WeightingResults1.mat';
filename_weight='WeightingResults_n0_66667_backup2.mat';
% filename_weight='WeightingResults_n2_5.mat';
load(filename_weight);

bay_resid=real(bay_resid);
lp_resid=real(lp_resid);
chi_sq_outletsfree=real(chi_sq_outletsfree);
v_s=real(v_s);
v_n=real(v_n);

sig_bay_weights=[-6.3 -6.4 -6.5];
sig_loma_weights=[-6.1 -6.4 -6.8];


bay_resid_weights=interp2(log_w_bay,log_w_lp,bay_resid,sig_bay_weights,sig_loma_weights,'nearest');
lp_resid_weights=interp2(log_w_bay,log_w_lp,lp_resid,sig_bay_weights,sig_loma_weights,'nearest');
chi_sq_outletsfree_weights=interp2(log_w_bay,log_w_lp,chi_sq_outletsfree,sig_bay_weights,sig_loma_weights,'nearest');
v_s_weights=interp2(log_w_bay,log_w_lp,v_s,sig_bay_weights,sig_loma_weights,'nearest');
v_n_weights=interp2(log_w_bay,log_w_lp,v_n,sig_bay_weights,sig_loma_weights,'nearest');


%%
figure
imagesc(log_w_bay,log_w_lp,chi_sq_outletsfree,'AlphaData',(~isnan(chi_sq_outletsfree))),
axis equal, axis xy, set(gca,'Ydir','normal')
% CT=cbrewer('seq', 'PuBuGn', 100); colormap(CT); colorbar
CT=cbrewer('div', 'Spectral', 100); colormap(CT); colorbar
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Reduced Chi-Squared Values Elevation Points')
font(14)
set(gcf,'Position',[10, 10, 600, 600])

figure
imagesc(log_w_bay,log_w_lp,abs(bay_resid.*1e3),'AlphaData',(~isnan(bay_resid))),
axis equal, axis xy, set(gca,'Ydir','normal')
% CT=cbrewer('seq', 'YlOrBr', 100); colormap(CT); colorbar
CT=cbrewer('div', 'BrBG', 100); colormap(CT); colorbar
% caxis([0 0.5])
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Bay Residuals [mm/yr]')
font(14)
set(gcf,'Position',[10, 10, 600, 600])

figure
imagesc(log_w_bay,log_w_lp,abs(lp_resid.*1e3),'AlphaData',(~isnan(bay_resid))),
axis equal, axis xy, set(gca,'Ydir','normal')
CT=cbrewer('seq', 'YlGn', 100); colormap(CT); colorbar
% CT=cbrewer('div', 'RdYlBu', 100); colormap(CT); colorbar
% caxis([0 0.5])
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Loma Prieta Residuals [mm/yr]')
font(14)
set(gcf,'Position',[10, 10, 600, 600])
%%
figure
imagesc(log_w_bay,log_w_lp,v_s.*1e3,'AlphaData',(~isnan(v_s))),
axis equal, axis xy, set(gca,'Ydir','normal')
CT=cbrewer('seq', 'Reds', 100); colormap(CT); colorbar
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Shear Plate Motion [mm/yr]')
font(14)
set(gcf,'Position',[10, 10, 600, 600])

figure
imagesc(log_w_bay,log_w_lp,v_n.*-1e3,'AlphaData',(~isnan(v_n))),
axis equal, axis xy, set(gca,'Ydir','normal')
caxis([-max(abs(v_n(:).*1e3)) max(abs(v_n(:).*1e3))])
CT=cbrewer('div', 'RdYlBu', 100); colormap(CT); colorbar
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Normal Plate Motion [mm/yr]')
font(14)
set(gcf,'Position',[10, 10, 600, 600])


figure
imagesc(log_w_bay,log_w_lp,k,'AlphaData',(~isnan(k))),
axis equal, axis xy, set(gca,'Ydir','normal')
% caxis([min(k(:)) max(k(:))])
CT=cbrewer('div', 'PiYG', 100); colormap(CT); colorbar
label('log \sigma Bay [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
hold on
plot(sig_bay_weights(1),sig_loma_weights(1),'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(2),sig_loma_weights(2),'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot(sig_bay_weights(3),sig_loma_weights(3),'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.05 max(log_w_bay)+.05])
ylim([min(log_w_lp)-.05 max(log_w_lp)+.05])
title('Log10(Erodibility) [m^{0.2}/yr]')
font(14)
set(gcf,'Position',[10, 10, 600, 600])

%%
figure
surf(log_w_bay,log_w_lp,chi_sq_outletsfree), shading flat
axis xy
CT=cbrewer('div', 'Spectral', 100); colormap(CT); colorbar
% label('log \sigma Loma Prieta [mm/yr]','log \sigma Loma Prieta [mm/yr]',14)
xlabel('log \sigma Bay [mm/yr]')
ylabel('log \sigma Loma Prieta [mm/yr]')
zlabel('Red Chi-Sq')
hold on
plot3(sig_bay_weights(1),sig_loma_weights(1),chi_sq_outletsfree_weights(1)+.2,'ok','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot3(sig_bay_weights(2),sig_loma_weights(2),chi_sq_outletsfree_weights(2)+.2,'dk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',12)
plot3(sig_bay_weights(3),sig_loma_weights(3),chi_sq_outletsfree_weights(3)+.2,'pk','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20)
hold off
xlim([min(log_w_bay)-.01 max(log_w_bay)+.01])
ylim([min(log_w_lp)-.01 max(log_w_lp)+.01])
title('Reduced Chi-Squared Values Elevation Points')
font(14)
set(gcf,'Position',[10, 10, 1100, 600])
view(125,50)

%%
%%

figure
pcolor(LOG_SIG_LOMA_SPACE,LOG_SIG_BAY_SPACE,LOG_WRSS_elev_atminWRSScombined_SPACE), shading flat
% CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
CT=cbrewer('seq', 'PuBu', 100); colormap(CT)
colorbar
% colormap(jet), colorbar
label('log \sigma Loma Prieta [mm/yr]','log \sigma bay [mm/yr]',14)
hold on
plot(sig_loma_weights(1),sig_bay_weights(3),'pk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',20)
plot(sig_loma_weights(2),sig_bay_weights(2),'dk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',12)
plot(sig_loma_weights(3),sig_bay_weights(1),'ok','LineWidth',2,'MarkerFaceColor','r','MarkerSize',12)
% title(['WRSS elevation, \sigma = ',num2str(sig_elev),' [m], K=10^{',num2str(log_k_model(1)),'}'])
title('Log WRSS Elevation at Min Combined WRSS')
axis equal
xlim([-3.9 -3.3])
ylim([-4.1 -3.3])
hold off
font(14)
set(gcf,'Position',[10, 10, 550, 600])

xlabel('Right-lateral shear [mm/yr]')
xlabel('Right-lateral shear [mm/yr]')
title('Elevation+Bay Shore+Loma Prieta')

%%
figure
imagesc(log_w_bay,log_w_lp,bay_resid), axis equal, set(gca,'Ydir','normal')


CT=cbrewer('div', 'Spectral', 50); colormap(CT)
    colorbar