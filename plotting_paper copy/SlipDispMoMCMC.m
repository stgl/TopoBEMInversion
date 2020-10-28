clear all
close all

% addpath(genpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries'))
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries')
addpath(genpath('/Users/felipearon/Documents/POSTDOC/Codes/BEM/Jacks_tribemx/tribem2018'))
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh')
addpath('/Users/felipearon/Documents/POSTDOC/Codes/MCMC/MCMC_postprocess')
addpath(genpath('/Users/felipearon/Documents/Codes/colorpalettes/'))

chi_file='obs_data_chi.txt';
grid_file='obs_data_grid.txt';


%in mm/yr
%weighting a: Bay= 10^{-3.95}, Loma Prieta= 10^{-3.65} (weighting a)
%weighting i: Bay= 10^{-3.75}, Loma Prieta= 10^{-3.45}
weight_flag='aaa';
% weight_flag='eee';
% weight_flag='iii';
ML_l95_u95_push=load([weight_flag,'_weight_ML_l95_u95_push.txt']);
ML_l95_u95_shear=load([weight_flag,'_weight_ML_l95_u95_shear.txt']);

%in mm/yr, 
MCMC_ML_shear_push=[ML_l95_u95_shear(1),ML_l95_u95_push(1)];
MCMC_l95_shear_push=[ML_l95_u95_shear(2),ML_l95_u95_push(2)];
MCMC_u95_shear_push=[ML_l95_u95_shear(3),ML_l95_u95_push(3)];

% Control the model solution (ML, u95 or l95)
% MCMC_shear_push=MCMC_ML_shear_push;
% MCMC_shear_push=MCMC_l95_shear_push;
MCMC_shear_push=MCMC_u95_shear_push;

% %in mm/yr, Bay= 10^{-3.95}, Loma Prieta= 10^{-3.65} (weighting a)
% MCMC_MaxLike_shear_push=[0.016164034616631,-0.000829645383525]*1e3;

% %in mm/yr, Bay= 10^{-3.95}, Loma Prieta= 10^{-3.45} (weighting c)
% MCMC_MaxLike_shear_push=[,]*1e3;

% %in mm/yr, Bay= 10^{-3.85}, Loma Prieta= 10^{-3.55} (weighting e)
% MCMC_MaxLike_shear_push=[0.012676276042513,0.000841995564879]*1e3;

% %in mm/yr, Bay= 10^{-3.75}, Loma Prieta= 10^{-3.65} (weighting g)
% MCMC_MaxLike_shear_push=[,]*1e3;

% %in mm/yr, Bay= 10^{-3.75}, Loma Prieta= 10^{-3.45} (weighting i)
% MCMC_MaxLike_shear_push=[0.010511571398565,0.001758685369214]*1e3;

cl='2.5';
% v=3;

xyz_chi=load(chi_file);
xyz_grid=load(grid_file);

len_chi=length(xyz_chi);
len_chi_str=num2str(len_chi);
len_grid=length(xyz_grid);
len_grid_str=num2str(len_grid);

xyz=[xyz_chi;xyz_grid];
    
% obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);

% Faults
BE_filename= ['BE_cl',cl,'.msh'];
SHAMTV_filename= ['SHAMTV_cl',cl,'.msh'];
SA_filename= ['SA_cl',cl,'.msh'];
NAdecoll_filename= ['NA_decoll_cl',cl,'.msh'];
PAdecoll_filename= ['PA_decoll_cl',cl,'.msh'];

% Box patches (tectonic BCs)
BoxL_filename= ['box_L_cl',cl,'.msh'];
BoxR_filename= ['box_R_cl',cl,'.msh'];
BoxLL_filename= ['box_LL_cl',cl,'.msh'];
BoxUR_filename= ['box_UR_cl',cl,'.msh'];
BoxUL_filename= ['box_UL_cl',cl,'.msh'];
BoxLR_filename= ['box_LR_cl',cl,'.msh'];

nr_box_patches=6;

faultnames=[BE_filename,' ',SHAMTV_filename,' ',SA_filename, ' ',...
    NAdecoll_filename, ' ',PAdecoll_filename,' ',BoxL_filename,' ',...
    BoxR_filename,' ',BoxLL_filename,' ',BoxUR_filename,' ',...
    BoxUL_filename,' ',BoxLR_filename];

faults=ReadPatches(faultnames); % in km

ends = cumsum(faults.nEl); % Give ending indices of the 3 faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the 3 faults

centroids=PatchCentroid(faults.c,faults.v); % xyz of faults centroids, km ENU

inname_shear=['GreenFncs_cl',cl,'_shear.mat']; load(inname_shear);
inname_push=['GreenFncs_cl',cl,'_push.mat']; load(inname_push);

u_shear=reshape(obs_shear.u,3,[])';
u_shear_chi=u_shear(1:len_chi,3);
u_shear_grid=u_shear(len_chi+1:end,3);

u_push=reshape(obs_push.u,3,[])';
u_push_chi=u_push(1:len_chi,3);
u_push_grid=u_push(len_chi+1:end,3);

%%
str_MCMC=[slip_shear(:,1) slip_push(:,1)]*(MCMC_shear_push'./2);
% str_MCMC=[slip_shear(:,1) slip_push(:,1)]*(MCMC_ML_shear_push');
dip_MCMC=[slip_shear(:,2) slip_push(:,2)]*(MCMC_shear_push'./2);
% dip_MCMC=[slip_shear(:,2) slip_push(:,2)]*(MCMC_ML_shear_push');
nor_MCMC=[slip_shear(:,3) slip_push(:,3)]*(MCMC_shear_push'./2);
% nor_MCMC=[slip_shear(:,3) slip_push(:,3)]*(MCMC_ML_shear_push');
slip_total_MCMC=sqrt((str_MCMC.^2)+(dip_MCMC.^2));


slip_xyz=slip2cart(faults.c.*1e3,faults.v,str_MCMC,dip_MCMC,nor_MCMC);
%%


ux_MCMC=[u_shear(:,1) u_push(:,1)]*(MCMC_shear_push'./2);
% ux_MCMC=[u_shear(:,1) u_push(:,1)]*(MCMC_ML_shear_push');
ux_MCMC_chi=ux_MCMC(1:len_chi);
ux_MCMC_grid=ux_MCMC(len_chi+1:end);

uy_MCMC=[u_shear(:,2) u_push(:,2)]*(MCMC_shear_push'./2);
% uy_MCMC=[u_shear(:,2) u_push(:,2)]*(MCMC_ML_shear_push');
uy_MCMC_chi=uy_MCMC(1:len_chi);
uy_MCMC_grid=uy_MCMC(len_chi+1:end);

uz_MCMC=[u_shear(:,3) u_push(:,3)]*(MCMC_shear_push'./2);
% uz_MCMC=[u_shear(:,3) u_push(:,3)]*(MCMC_ML_shear_push');
uz_MCMC_chi=uz_MCMC(1:len_chi);
uz_MCMC_grid=uz_MCMC(len_chi+1:end);

% bay_pts=[582,4147;592,4137];


%%


%in m/yr, Pa; out Nm/yr, m^2
ShearMod=3.2e10;
[momentrate,elemarea]=MomentRate(faults.c.*1e3,faults.v,slip_total_MCMC.*1e-3,ShearMod);

mo_rate_berrocal=sum(momentrate(begs(1):ends(1)));
mo_rate_shannon=sum(momentrate(begs(2):ends(2)));
mo_rate_FTB=mo_rate_berrocal+mo_rate_shannon;

mw_rate_berrocal=(2/3)*(log10(mo_rate_berrocal*1e7))-10.7;
mw_rate_shannon=(2/3)*(log10(mo_rate_shannon*1e7))-10.7;
mw_rate_FTB=(2/3)*(log10(mo_rate_FTB*1e7))-10.7;

%%
area_aver_slip_berrocal=sum((slip_total_MCMC(begs(1):ends(1))).*(elemarea(begs(1):ends(1))./sum(elemarea(begs(1):ends(1)))));
area_aver_slip_shannon=sum((slip_total_MCMC(begs(2):ends(2))).*(elemarea(begs(2):ends(2))./sum(elemarea(begs(2):ends(2)))));
max_slip_berrocal=max(slip_total_MCMC(begs(1):ends(1)));
max_slip_shannon=max(slip_total_MCMC(begs(2):ends(2)));
area_aver_slip_berrocal_shannon=sum((slip_total_MCMC(begs(1):ends(2))).*(elemarea(begs(1):ends(2))./sum(elemarea(begs(1):ends(2)))));

%%
load('BE_trace.txt');
load('SHAMTV_trace.txt');
load('SA_trace.txt');
load('box_trace.txt');

nan_bk=[NaN,NaN];
traces=[BE_trace;nan_bk;SHAMTV_trace;nan_bk;SA_trace;nan_bk;box_trace];

load('study_area_paper_utm.txt');
study_area_paper_utm=study_area_paper_utm.*1e-3;
%%
fault_colors=ones(ends(5),1);
fault_colors(1:ends(1))=fault_colors(1:ends(1)).*1;
fault_colors(begs(2):ends(2))=fault_colors(begs(2):ends(2)).*2;
fault_colors(begs(3):ends(3))=fault_colors(begs(3):ends(3)).*3;
fault_colors(begs(4):ends(5))=fault_colors(begs(4):ends(5)).*4;


%% Paper ready (for Fig.1)
% meshview(faults.c, faults.v(begs(1):ends(5),:));
meshview(faults.c, faults.v(begs(1):ends(5),:),fault_colors);
% colorbar, colormap(redblue), caxis([min(str_MCMC(begs(3):ends(3))), min(str_MCMC(begs(3):ends(3)))*-1])
CT=cbrewer('qual', 'Set3', 4); colormap(flipud(CT))
axis equal
% hold on
% plot3(study_area_paper_utm(:,1),study_area_paper_utm(:,2),study_area_paper_utm(:,3),'-k','LineWidth',0.5)
view(52,9)
% shading flat
xlim([384,837.5])
ylim([3883,4350])
zlim([-23 0])

zoom(9);

xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
% box on
% grid on
title('Geometric Arrangment')
% hold off
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
% set(gcf,'Position',[0, 0, 1200, 180])
set(gcf,'Position',[0, 0, 1400, 600])

%% Paper ready Fig S1A
% meshview(faults.c, faults.v(begs(1):ends(5),:));
meshview(faults.c, faults.v(begs(1):ends(5),:),fault_colors);
% colorbar, colormap(redblue), caxis([min(str_MCMC(begs(3):ends(3))), min(str_MCMC(begs(3):ends(3)))*-1])
CT=cbrewer('qual', 'Set3', 4); colormap(flipud(CT))
axis equal
% hold on
% plot3(study_area_paper_utm(:,1),study_area_paper_utm(:,2),study_area_paper_utm(:,3),'-k','LineWidth',0.5)
view(0,20)
% shading flat
xlim([384,837.5])
ylim([3883,4350])
zlim([-23 0])

% box on
% grid on

zoom(2);

xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
% box on
% grid on
title('Geometric Arrangment')
% hold off
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
% set(gcf,'Position',[0, 0, 1200, 180])
set(gcf,'Position',[0, 0, 1300, 700])


%% Paper ready Fig S1B
% meshview(faults.c, faults.v(begs(1):ends(5),:));
meshview(faults.c, faults.v(begs(1):ends(5),:),fault_colors);
% colorbar, colormap(redblue), caxis([min(str_MCMC(begs(3):ends(3))), min(str_MCMC(begs(3):ends(3)))*-1])
CT=cbrewer('qual', 'Set3', 4); colormap(flipud(CT))
axis equal
% hold on
% plot3(study_area_paper_utm(:,1),study_area_paper_utm(:,2),study_area_paper_utm(:,3),'-k','LineWidth',0.5)
view(55,8)
% shading flat
xlim([384,837.5])
ylim([3883,4350])
zlim([-23 0])

% box on
% grid on

zoom(15);

xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
% box on
% grid on
title('Geometric Arrangment')
% hold off
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
% set(gcf,'Position',[0, 0, 1200, 180])
set(gcf,'Position',[0, 0, 800, 700])

%%
meshview(faults.c, faults.v(begs(1):ends(5),:), dip_MCMC(begs(1):ends(5)));
colorbar, colormap(jet), caxis([-2 2]), axis equal

meshview(faults.c, faults.v(begs(1):ends(5),:), str_MCMC(begs(1):ends(5)));
colorbar, colormap(jet)
caxis([-20 20])
axis equal

%%

meshview(faults.c, faults.v(begs(1):ends(3),:), slip_total_MCMC(begs(1):ends(3)));
colorbar, colormap(cool)
caxis([0 30])
shading flat
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Total-slip on San Andreas Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

meshview(faults.c, faults.v(begs(1):ends(3),:), str_MCMC(begs(1):ends(3)));
colorbar, colormap(jet)
caxis([-30 30])
shading flat
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Strike-slip on San Andreas Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

%%

%%% Berrocal
meshview(faults.c, faults.v(begs(1):ends(1),:), dip_MCMC(begs(1):ends(1)));
colorbar, colormap(jet)
shading flat
caxis([-1 1])
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Dip-slip on Berrocal Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

meshview(faults.c, faults.v(begs(1):ends(1),:), str_MCMC(begs(1):ends(1)));
colorbar, colormap(jet)
shading flat
caxis([-2 2])
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Strike-slip on Berrocal Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

meshview(faults.c, faults.v(begs(1):ends(1),:), slip_total_MCMC(begs(1):ends(1)));
% colorbar, colormap(flipud(hot))
colorbar, colormap(jet)

shading flat
caxis([0 2])

axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Total-slip on Berrocal Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

%% Berrocal paper ready
meshview(faults.c, faults.v(begs(1):ends(1),:), slip_total_MCMC(begs(1):ends(1)));
% colorbar, colormap(flipud(hot))
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar
% colorbar, cptcmap('GMT_no_green','flip',false,'ncol',100)
shading flat
caxis([0 2])

hold on
hq=quiver3(centroids(begs(1):ends(1),1),centroids(begs(1):ends(1),2),centroids(begs(1):ends(1),3),...
    slip_xyz(begs(1):ends(1),1),slip_xyz(begs(1):ends(1),2),slip_xyz(begs(1):ends(1),3),0.25,'ok',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','w',...
    'MarkerFaceColor','w',...
    'MarkerSize',3);

axis equal
% xlim([545,660])
% ylim([4050,4160])
% zlim([-21 0])
xlim([565,650])
ylim([4065,4145])
zlim([-15 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Total-slip on Berrocal Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-45 6])
% view([-44 10])
% set(gcf,'Position',[0, 0, 1100, 300])
set(gcf,'Position',[0, 400, 1400, 300])

%%

%%%%%%%
% Change arrow style quiver 2D using arrow annotation. Does not work in 3D
% Based on
% https://stackoverflow.com/questions/18776172/in-matlab-how-do-i-change-the-arrow-head-style-in-quiver-plot
%%%%%%


%%
%%% Shannon
meshview(faults.c, faults.v(begs(2):ends(2),:), dip_MCMC(begs(2):ends(2)));
colorbar, colormap(jet)
shading flat
caxis([-2 2])
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Dip-slip on Shannon-Monte Vista Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])


meshview(faults.c, faults.v(begs(2):ends(2),:), str_MCMC(begs(2):ends(2)));
colorbar, colormap(jet)
shading flat
caxis([-3 3])
axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]') 
zlabel('Depth [km]')
title('Modeled Strike-slip on Shannon-Monte Vista Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

meshview(faults.c, faults.v(begs(2):ends(2),:), slip_total_MCMC(begs(2):ends(2)));
% colorbar, colormap(flipud(hot))
colorbar, colormap(jet)
shading flat
caxis([0 2.5])

axis equal
xlim([545,660])
ylim([4050,4160])
zlim([-21 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Total-slip on Shannon-Monte Vista Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-44 10])
set(gcf,'Position',[0, 0, 1100, 300])

%% Shannon paper ready
meshview(faults.c, faults.v(begs(2):ends(2),:), slip_total_MCMC(begs(2):ends(2)));
% colorbar, colormap(flipud(hot))
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar
shading flat
caxis([0 2.5])

hold on
hq=quiver3(centroids(begs(2):ends(2),1),centroids(begs(2):ends(2),2),centroids(begs(2):ends(2),3),...
    slip_xyz(begs(2):ends(2),1),slip_xyz(begs(2):ends(2),2),slip_xyz(begs(2):ends(2),3),0.35,'ok',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','w',...
    'MarkerFaceColor','w',...
    'MarkerSize',3);

axis equal
% xlim([545,660])
% ylim([4050,4160])
% zlim([-21 0])
xlim([565,650])
ylim([4065,4145])
zlim([-15 0])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Modeled Total-slip on Shannon-Monte Vista Fault [mm/yr]')
font(16)
box on
grid on
hold off
view([-45 6])
% view([-44 10])
% set(gcf,'Position',[0, 0, 1100, 300])
set(gcf,'Position',[0, 0, 1400, 300])

%%
figure
% quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear(len_chi+1:end,1),u_shear(len_chi+1:end,2),1,'-mo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
quiver(xyz_grid(:,1),xyz_grid(:,2),ux_MCMC_grid,uy_MCMC_grid,1,'m')
hold on
plot(xyz_grid(:,1),xyz_grid(:,2),'og','MarkerFaceColor','k','MarkerSize',1.5)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid xy-disp MCMC MaxLike finest model (18,621 elements)')

%%
figure

% quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear(len_chi+1:end,1),u_shear(len_chi+1:end,2),1,'-mo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
scatter(xyz_grid(:,1),xyz_grid(:,2),50,uz_MCMC_grid,'filled')
colorbar, colormap(jet)
caxis([-1, 1])
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid uz-disp MCMC MaxLike finest model (18,621 elements)')

%%
figure

scatter(xyz_chi(:,1),xyz_chi(:,2),3,uz_MCMC_chi,'filled')
hold on
plot(traces(:,1),traces(:,2),'-k')
colorbar, colormap(jet), caxis([0 max(uz_MCMC_chi)]), axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')

%% Uplift rate paper ready
figure

scatter(xyz_chi(:,1),xyz_chi(:,2),8,uz_MCMC_chi,'filled')
% hold on
% plot(traces(:,1),traces(:,2),'-k')
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar

% caxis([0 max(uz_MCMC_chi)]), axis equal
caxis([0 1]), axis equal
xlim([560,640])
ylim([4080,4150])
% hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
font(16)
box on
% hold off
set(gcf,'Position',[0, 0, 900, 750])



%%

figure

scatter(xyz_chi(:,1),xyz_chi(:,2),3,uz_MCMC_chi,'filled')
hold on
plot(traces(:,1),traces(:,2),'-k')
colorbar, colormap(jet), caxis([0 max(uz_MCMC_chi)]), axis equal
xlim([550,640])
ylim([4075,4150])
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Modeled Uz (uplift-rate) [mm/yr]')
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])



%%

meshview(faults.c, faults.v(begs(1):ends(5),:), slip_shear(begs(1):ends(5), 1)); colorbar, colormap(redblue), axis equal
meshview(faults.c, faults.v(begs(1):ends(5),:), slip_shear(begs(1):ends(5), 2)); colorbar, colormap(redblue), axis equal

meshview(faults.c, faults.v(begs(1):ends(5),:), slip_push(begs(1):ends(5), 1)); colorbar, colormap(redblue), axis equal
meshview(faults.c, faults.v(begs(1):ends(5),:), slip_push(begs(1):ends(5), 2)); colorbar, colormap(redblue), axis equal

%%
meshview(faults.c, faults.v(begs(3):ends(3),:), slip_shear(begs(3):ends(3), 1)); colorbar, colormap(redblue), axis equal
meshview(faults.c, faults.v(begs(3):ends(3),:), slip_shear(begs(3):ends(3), 2)); colorbar, colormap(redblue), axis equal

%%
figure

subplot(1,2,1)
% quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear(len_chi+1:end,1),u_shear(len_chi+1:end,2),1,'-mo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear_grid(:,1),u_shear_grid(:,2),1,'m')
hold on
plot(xyz_grid(:,1),xyz_grid(:,2),'og','MarkerFaceColor','k','MarkerSize',1.5)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid xy-disp Shear Green Fcns finest model (18,621 elements)')

subplot(1,2,2)
quiver(xyz_grid(:,1),xyz_grid(:,2),u_push_grid(:,1),u_push_grid(:,2),1,'m')
hold on
plot(xyz_grid(:,1),xyz_grid(:,2),'og','MarkerFaceColor','k','MarkerSize',1.5)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid xy-disp Push Green Fcns finest model (18,621 elements)')

%%
figure

subplot(1,2,1)
% quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear(len_chi+1:end,1),u_shear(len_chi+1:end,2),1,'-mo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
scatter(xyz_grid(:,1),xyz_grid(:,2),50,u_shear_grid(:,3),'filled')
colorbar, colormap(jet), caxis([0, 0.18])
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid uz-disp Shear Green Fcns finest model (18,621 elements)')

subplot(1,2,2)
scatter(xyz_grid(:,1),xyz_grid(:,2),50,u_push_grid(:,3),'filled')
colorbar, colormap(jet), caxis([0, 0.18])
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Grid uz-disp Push Green Fcns finest model (18,621 elements)')


%%
% figure, scatter(xyz_grid(:,1),xyz_grid(:,2),10,u_shear(len_chi+1:end,3),'filled'), colorbar, colormap(jet), caxis([0 0.2]), axis equal
figure

subplot(1,2,1)
scatter(xyz_chi(:,1),xyz_chi(:,2),3,u_shear(1:len_chi,3),'filled')
hold on
plot(traces(:,1),traces(:,2),'-k')
colorbar, colormap(jet), caxis([0 max(u_shear(1:len_chi,3))]), axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Uz (uplift) chi obs points Shear Green Fcns finest model (18,621 elements)')


% figure, scatter(xyz_grid(:,1),xyz_grid(:,2),10,u_push(len_chi+1:end,3),'filled'), colorbar, colormap(jet), caxis([0 0.2]), axis equal
subplot(1,2,2)
scatter(xyz_chi(:,1),xyz_chi(:,2),3,u_push(1:len_chi,3),'filled')
hold on
plot(traces(:,1),traces(:,2),'-k')
colorbar, colormap(jet), caxis([0 max(u_push(1:len_chi,3))]), axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Uz (uplift) chi obs points Push Green Fcns finest model (18,621 elements)')


%%

G_tect=[u_shear(1:len_chi,3),u_push(1:len_chi,3)];
dlmwrite('G_tect.txt',G_tect);