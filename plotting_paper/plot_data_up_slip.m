%% Input inversion results

scenario = 'weak'; % 'weak', 'medium', 'strong'

if strcmp(scenario,'strong')
    title_analysis='"Strong" Weighting Model';
elseif strcmp(scenario,'medium')
    title_analysis='"Medium" Weighting Model';
else
    title_analysis='"Weak" Weighting Model';
end

simple_model=0;


% % Setting the correct directories
% p=pathdef; path(p)


addpath(genpath('/Users/felipearon/Documents/Codes/colorpalettes'))
addpath(genpath('/Users/felipearon/Documents/Codes/utils'))
addpath(genpath('/Users/felipearon/Documents/GitHub/TopoBEMInversion'))

% addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries')
addpath(genpath('/Users/felipearon/Documents/POSTDOC/Codes/BEM/Jacks_tribemx/tribem2018'))
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh')

%Check
% addpath('/Users/felipearon/Documents/POSTDOC/Codes/MCMC/MCMC_postprocess')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation grids
chi_file='obs_data_chi.txt';
grid_file='obs_data_grid.txt';

xyz_chi=load(chi_file);
xyz_grid=load(grid_file);

len_chi=length(xyz_chi);
len_chi_str=num2str(len_chi);
len_grid=length(xyz_grid);
len_grid_str=num2str(len_grid);

xyz=[xyz_chi;xyz_grid];
 
% obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);

load('study_area_paper_utm.txt');
study_area_paper_utm=study_area_paper_utm.*1e-3;

% Fault traces
load('BE_trace.txt');
load('SHAMTV_trace.txt');
load('SA_trace.txt');
load('box_trace.txt');

nan_bk=[NaN,NaN];
traces=[BE_trace;nan_bk;SHAMTV_trace;nan_bk;SA_trace;nan_bk;box_trace];


%% Load inversion results

if simple_model
    
    load(strcat('results/singleK_', scenario, '_enforcement'));
    load(strcat('input_files/m_', scenario, '_simple.txt'));
    load(strcat('input_files/m_', scenario, '_std_simple.txt'));
    
    m_file=strcat('m_',scenario,'_simple.txt'); m_ML=load(m_file);
    m_std_file=strcat('m_',scenario,'_std_simple.txt'); m_std=load(m_std_file);
    
else
    
    load(strcat('results/multiK_', scenario, '_enforcement'));
    load(strcat('input_files/m_', scenario, '.txt'));
    load(strcat('input_files/m_', scenario, '_std.txt'));
    
    m_file=strcat('m_',scenario,'.txt'); m_ML=load(m_file);
    m_std_file=strcat('m_',scenario,'_std.txt'); m_std=load(m_std_file);
    
    load('geo_map.txt');
    nr_nounits=8;
    
    K_nr=NaN(n_K,1);
    for i=1:n_K
        if i==1
            K_nr(i,1)=geo_map(i,3);
        else
            K_nr(i,1)=geo_map(i,3)-geo_map(i-1,3);
        end
        
        
    end
    
    keys = textread('geo_keys.txt','%s');
    fid = fopen('geokeys.txt');
    geocol = textscan(fid,'%s%d%d%f%f%f%s%d');
    fclose(fid);
    
    plotorder = double(geocol{8}); 
    
    % K values geo scale ordered
    
    K_geo=K(plotorder)';
    K_geo_std=m_std(3:end);K_geo_std=K_geo_std(plotorder);
    
    K_geo_nodummy=K_geo(1:end-nr_nounits);
    K_geo_nodummy_std=K_geo_std(1:end-nr_nounits);
    K_geo_nodummy_2std=K_geo_nodummy_std.*2;
    
    K_geo_nr=K_nr(plotorder);
    K_geo_nodummy_nr=K_geo_nr(1:end-nr_nounits);
    
    % Geologic units colors for elevation points
    col=ones(nr_data,3)./255;
    % col(geo_map(i,2):geo_map(i,3),3).*geocol{6}(i)]./255;
    
    for i=1:n_K-nr_nounits
        col(geo_map(i,2):geo_map(i,3),:) =...
            [col(geo_map(i,2):geo_map(i,3),1).*geocol{4}(i) ...
            col(geo_map(i,2):geo_map(i,3),2).*geocol{5}(i) ...
            col(geo_map(i,2):geo_map(i,3),3).*geocol{6}(i)];
    end
    
end


ML_shear_push=[vshear, vconverge]; %in m
STD_shear_push=[m_std(1)/1e3, m_std(2)/1e3]; %in m
STD2_shear_push=2*STD_shear_push; %in m

%% Load BEM model Green's functions of surface deformation and fault slip

cl='2.5'; % Characteristic length mesh

% Faults meshes simple and complex model
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

if simple_model
    run_flag='a';
    Date_PID_flag='2019_4_13_PID1880';
    
    % Faults meshes FTB
    BE_filename= ['BE_simple_cl',cl,'.msh'];
    SHAMTV_filename= ['SHAMTV_simple_cl',cl,'.msh'];
    
    inname_shear=[run_flag,'_output_SCM_Greens_simple_shear_',Date_PID_flag]; load(inname_shear);
    inname_push=[run_flag,'_output_SCM_Greens_simple_push_',Date_PID_flag]; load(inname_push);
    
else
    
    % Faults meshes FTB
    BE_filename= ['BE_cl',cl,'.msh'];
    SHAMTV_filename= ['SHAMTV_cl',cl,'.msh'];
    
    inname_shear=['GreenFncs_cl',cl,'_shear.mat']; load(inname_shear);
    inname_push=['GreenFncs_cl',cl,'_push.mat']; load(inname_push);
    
end




faultnames=[BE_filename,' ',SHAMTV_filename,' ',SA_filename, ' ',...
    NAdecoll_filename, ' ',PAdecoll_filename,' ',BoxL_filename,' ',...
    BoxR_filename,' ',BoxLL_filename,' ',BoxUR_filename,' ',...
    BoxUL_filename,' ',BoxLR_filename];

faults=ReadPatches(faultnames); % in km

ends = cumsum(faults.nEl); % Give ending indices of the 3 faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the 3 faults

centroids=PatchCentroid(faults.c,faults.v); % xyz of faults centroids, km ENU


u_shear=reshape(obs_shear.u,3,[])';
u_shear_chi=u_shear(1:len_chi,3);
u_shear_grid=u_shear(len_chi+1:end,3);

u_push=reshape(obs_push.u,3,[])';
u_push_chi=u_push(1:len_chi,3);
u_push_grid=u_push(len_chi+1:end,3);


%%  Compute maximum likelihood surface displacement and fault slip

% Surface displacement rates all in [m]
ux_model=[u_shear(:,1) u_push(:,1)]*(ML_shear_push'./2);
ux_model_chi=ux_model(1:len_chi);
ux_model_grid=ux_model(len_chi+1:end);

uy_model=[u_shear(:,2) u_push(:,2)]*(ML_shear_push'./2);
uy_model_chi=uy_model(1:len_chi);
uy_model_grid=uy_model(len_chi+1:end);

uz_model=[u_shear(:,3) u_push(:,3)]*(ML_shear_push'./2);
uz_model_chi=uz_model(1:len_chi);
uz_model_grid=uz_model(len_chi+1:end);

% bay_pts=[582,4147;592,4137];

% All units in [m]
mean_upliftrate=mean(uz_model_chi);


% Slip rates all in [m]
str_model=[slip_shear(:,1) slip_push(:,1)]*(ML_shear_push'./2);
% str_MCMC=[slip_shear(:,1) slip_push(:,1)]*(MCMC_ML_shear_push');
dip_model=[slip_shear(:,2) slip_push(:,2)]*(ML_shear_push'./2);
% dip_MCMC=[slip_shear(:,2) slip_push(:,2)]*(MCMC_ML_shear_push');
nor_model=[slip_shear(:,3) slip_push(:,3)]*(ML_shear_push'./2);
% nor_MCMC=[slip_shear(:,3) slip_push(:,3)]*(MCMC_ML_shear_push');
slip_total_model=sqrt((str_model.^2)+(dip_model.^2));


slip_xyz=slip2cart(faults.c.*1e3,faults.v,str_model,dip_model,nor_model);


%% Compute moment accrual rates

%in m/yr, Pa; out Nm/yr, m^2
ShearMod=3.2e10; % [Pa]
[momentrate,elemarea]=MomentRate(faults.c.*1e3,faults.v,slip_total_model,ShearMod);

mo_rate_berrocal=sum(momentrate(begs(1):ends(1)));
mo_rate_shannon=sum(momentrate(begs(2):ends(2)));
mo_rate_FTB=mo_rate_berrocal+mo_rate_shannon;

mw_rate_berrocal=(2/3)*(log10(mo_rate_berrocal*1e7))-10.7;
mw_rate_shannon=(2/3)*(log10(mo_rate_shannon*1e7))-10.7;
mw_rate_FTB=(2/3)*(log10(mo_rate_FTB*1e7))-10.7;

%
% All units in [m]
area_aver_slip_berrocal=sum((slip_total_model(begs(1):ends(1))).*(elemarea(begs(1):ends(1))./sum(elemarea(begs(1):ends(1)))));
area_aver_slip_shannon=sum((slip_total_model(begs(2):ends(2))).*(elemarea(begs(2):ends(2))./sum(elemarea(begs(2):ends(2)))));
max_slip_berrocal=max(slip_total_model(begs(1):ends(1)));
max_slip_shannon=max(slip_total_model(begs(2):ends(2)));
area_aver_slip_berrocal_shannon=sum((slip_total_model(begs(1):ends(2))).*(elemarea(begs(1):ends(2))./sum(elemarea(begs(1):ends(2)))));


%%
fault_colors=ones(ends(5),1);
fault_colors(1:ends(1))=fault_colors(1:ends(1)).*1;
fault_colors(begs(2):ends(2))=fault_colors(begs(2):ends(2)).*2;
fault_colors(begs(3):ends(3))=fault_colors(begs(3):ends(3)).*3;
fault_colors(begs(4):ends(5))=fault_colors(begs(4):ends(5)).*4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Paper ready Fig 2 inset
figure
histogram(residuals_elev_outletsfree./1e3,50,'FaceColor','k','EdgeColor',[.8 .8 .8])
% histogram(resi_elev_ML,100,'EdgeAlpha',0.5)

title([title_analysis, '\newline Hist Residuals Channel Elevation (Data-ML)'])
xlim([-0.5 0.5])
xlabel('Elevation [km]')
font(16)
set(gcf,'Position',[0, 0, 500, 700])


figure
histogram(outlets_resid./1e3,'FaceColor','y','EdgeColor',[.8 .8 .8])
% histogram(resi_elev_ML,100,'EdgeAlpha',0.5)

title([title_analysis, '\newline Hist Residuals Outlet Elevation (Data-ML)'])
xlim([-0.5 0.5])
xlabel('Elevation [km]')
font(16)
set(gcf,'Position',[0, 0, 500, 700])

%%
figure

subplot(1,2,1)
dscatter(e_chan,elev_pred_outletsfree)
% colorbar
colormap(jet)
axis equal
xlim([0 1100])
ylim([0 1100])
label('Channel elevations data','Predicted elevations Max Like')

subplot(1,2,2)
dscatter(log10(e_chan),log10(elev_pred_outletsfree))
% colorbar
colormap(jet)
axis equal
xlim([1.6 3.2])
ylim([1.6 3.2])
label('Log Channel elevations data','Log Predicted elevations Max Like')

set(gcf,'position',[150   372   750   350])
supertitle(['Dscatter plot Predicted vs Data ', title_analysis])
font(16)

set(gcf,'Position',[0, 0, 1200, 500])


%% Paper ready Fig 2
%2A
ax=figure;
scatter(x_channels./1e3,y_channels./1e3,8,e_chan,'filled')
caxis([0 1000])
% colorbar, colormap(parula)
colorbar, cptcmap('GMT_drywet',ax,'flip',true)
% hold on
% plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title(['Channel Elevations Data [m]'])
font(16)
box on
% hold off
set(gcf,'Position',[0, 0, 900, 750])

%Fig 2B
ax=figure;
scatter(x_channels./1e3,y_channels./1e3,10,elev_pred_outletsfree,'filled')
caxis([0 1000])
% colorbar, colormap(jet)
colorbar, cptcmap('GMT_drywet',ax,'flip',true)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title([title_analysis,', Channel Elev Maximum Likelihood [m]'])
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])


%Fig 2C
ax=figure;
scatter(x_channels./1e3,y_channels./1e3,8,residuals_elev_outletsfree,'filled')
caxis([-250 250])
% colorbar, colormap(ax,redblue)
% colorbar, cptcmap('GMT_polar',ax,'flip',false)
CT=cbrewer('div', 'RdGy', 100); colorbar, colormap(ax,CT)

hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title([title_analysis,', Residuals Channel Elevation (Data-ML) [m]'])
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])


%% Fig 3A Berrocal paper ready
meshview(faults.c, faults.v(begs(1):ends(1),:), slip_total_model(begs(1):ends(1)).*1e3);
% colorbar, colormap(flipud(hot))
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar
% colorbar, cptcmap('GMT_no_green','flip',false,'ncol',100)
shading flat
% caxis([0 2]) %strong
% caxis([0 1]) %weak
caxis([0 1.5]) %medium

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
title([title_analysis,', ML Total-slip on Berrocal Fault [mm/yr]'])
font(16)
box on
grid on
hold off
view([-45 6])
% view([-44 10])
% set(gcf,'Position',[0, 0, 1100, 300])
set(gcf,'Position',[0, 400, 1400, 300])

%% Shannon paper ready
meshview(faults.c, faults.v(begs(2):ends(2),:), slip_total_model(begs(2):ends(2)).*1e3);
% colorbar, colormap(flipud(hot))
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar
shading flat
% caxis([0 2.5]) %strong
% caxis([0 1.5]) %weak
caxis([0 2]) %medium

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
title([title_analysis,', ML Total-slip on Shannon-Monte Vista Fault [mm/yr]'])
font(16)
box on
grid on
hold off
view([-45 6])
% view([-44 10])
% set(gcf,'Position',[0, 0, 1100, 300])
set(gcf,'Position',[0, 0, 1400, 300])


%% Uplift rate paper ready
figure

scatter(xyz_chi(:,1),xyz_chi(:,2),8,uz_model_chi.*1e3,'filled')
hold on
plot(traces(:,1),traces(:,2),'-k')
% colorbar, colormap(jet)
CT=cbrewer('div', 'RdYlBu', 100); colormap(flipud(CT))
colorbar

% caxis([0 max(uz_MCMC_chi)]), axis equal
caxis([0 1])
axis equal
xlim([560,640])
ylim([4080,4150])
% hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title([title_analysis,', ML Uz (uplift rate) at Channel Observation Points'])
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])




%%
% Here down is recycle
%%    
figure
for k = 1:n_K
    i = find(ismember(plotorder,k));
    subplot(6,5,k); hold on
    col = [geocol{4}(i) geocol{5}(i) geocol{6}(i)]./255;
    h=area(Kxx(i,:),Kf(i,:),'facecolor',col,'linewidth',2);
%     area(Kxx(k,:),Kf(k,:),'facecolor','b','linewidth',2)
%         title(['log K_{' num2str(i) '}'])
    title([keys{i},' ',num2str((geo_map(i,3))-(geo_map(i,2))+1)])    %this one needs to be uncomment as well apr13/2018
%     title([keys(i),' ',num2str((geo_map(i,3))-(geo_map(i,2))+1)])    %this one needs to be uncomment as well apr13/2018
    font(16)
    
%         vline(Kml(i),'-r',2);
        vline(Kquant(i,1),'--r',2);
        vline(Kquant(i,2),'--r',2);
    
    vline(K_ML(i),'-r',2);
    
    
    if k == 1
        label('log K', 'Prob. density',16)
%     elseif k==2
%         title(title_analysis)
    end
    axis tight
%     xlim([-6 -1])
%     xlim([-5.5 -2.5])
    x_ext=abs((min(Kxx(i,:))-max(Kxx(i,:)))/2);
    xlim([min(Kxx(i,:))-x_ext/2 max(Kxx(i,:))+x_ext/2])
    xData = get(h,'XData');
    set(gca,'Xtick',linspace(xData(1)-x_ext/4,xData(end)+x_ext/4,3))
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.3f'))
    
    box on
    
end
hold off
% title(ax,title_analysis)
supertitle([title_analysis, '\newline '])
% supertitle('your title \newline ') 
font(16)

set(gcf,'Position',[0, 0, 1200, 850])


%% Paper ready Plot K BCs posteriors meaningful data

keys = textread('geo_keys.txt','%s');
fid = fopen('geokeys.txt');
geocol = textscan(fid,'%s%d%d%f%f%f%s%d');
fclose(fid);

% Kquant = Kstats.Kquant;
% Kxx = Kstats.Kxx;
% Kf = Kstats.Kf;
% Kml = Kstats.Kml;

plotorder = double(geocol{8});

nr_nounits=8;
    
figure
for k = 1:Ngeo-nr_nounits
    i = find(ismember(plotorder,k));
%     subplot(6,4,k); hold on
    subplot(7,3,k); hold on
    col = [geocol{4}(i) geocol{5}(i) geocol{6}(i)]./255;
    h=area(Kxx(i,:),Kf(i,:),'facecolor',col,'linewidth',2);
%     area(Kxx(k,:),Kf(k,:),'facecolor','b','linewidth',2)
%         title(['log K_{' num2str(i) '}'])
    title([keys{i},' ',num2str((geo_map(i,3))-(geo_map(i,2))+1)])    %this one needs to be uncomment as well apr13/2018
%     title([keys(i),' ',num2str((geo_map(i,3))-(geo_map(i,2))+1)])    %this one needs to be uncomment as well apr13/2018
    font(16)
    
%         vline(Kml(i),'-r',2);
        vline(Kquant(i,1),'--r',2);
        vline(Kquant(i,2),'--r',2);
    
    vline(K_ML(i),'-r',2);
    
    
    if k == 19
%     if k == 19 || k == 20 || k == 21
        label('log10 K [m^{0.2}/yr]', 'Prob. density',16)
        % xlabel('log K [m^{0.2}/yr]','FontSize',16)
  % elseif k == 1 || k == 4 || k == 7 || k == 10 || k == 13 || k == 16 || k == 19
        % ylabel('Prob. density','FontSize',16)
%     elseif k==2
%         title(title_analysis)
    end
    axis tight
%     xlim([-6 -1])
%     xlim([-5.5 -2.5])
    x_ext=abs((min(Kxx(i,:))-max(Kxx(i,:)))/2);
    xlim([min(Kxx(i,:))-x_ext/2 max(Kxx(i,:))+x_ext/2])    
    
    
    xData = get(h,'XData');
    set(gca,'Xtick',linspace(xData(1)-x_ext/4,xData(end)+x_ext/4,3))
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.3f'))
    
    box on
    
end
hold off
% title(ax,title_analysis)
supertitle({ '',title_analysis, '', '', ''})
% supertitle('your title \newline ') 
font(16)

set(gcf,'Position',[0, 0, 650, 800])


%% Save 

filename_model_logK=[weight_flag,'_weight_ML_l95_u95_logK_model.txt'];


fileID = fopen(filename_model_logK,'w');
fprintf(fileID,'%7.4f %7.4f %7.4f\n',[K_ML, Kstats.Kquant]');
fclose(fileID);

filename_model_logK=[weight_flag,'_weight_ML_l95_u95_logK_geo.txt'];
fileID = fopen(filename_model_logK,'w');
for k = 1:Ngeo
    i = find(ismember(plotorder,k));
    fprintf(fileID,'%7.4f %7.4f %7.4f\n',[K_ML(i), Kquant(i,1),Kquant(i,2)]');  
end
fclose(fileID);

%% Prepare inputs for topography obj function
theta = 0.4;
% Modified by Felipe May 24, 2017
% sig_elev=16;
% sig_elev=105.1636; %aaa
% sig_elev=84.3588; %eee
sig_elev=73.0495; %iii

[geo_map, d, Ginv_elev, channel_indexes, channel_elevations, n_K, compare_elements, litho_chan] = prepare_inputs_for_objective_function(theta);
channel_elevations_to_compare = channel_elevations(compare_elements);
%%
% Modified by Felipe March 8, 2018

G_tect=load('G_tect_shear_push.txt');
G_bay=load('G_bay_shear_push.txt');
G_lomaprieta=load('G_lomaprieta_shear_push.txt');

% Take account of the 2 BC due to 1 imposed on each side of the box
G_tect=G_tect./2; G_bay=G_bay./2; G_lomaprieta=G_lomaprieta./2;


v=0.25;
% v=0.49;

% % Far field slip G
% holocene_slip_SAfar=24;  %mm/yr, from Kaj's paper 

nr_bay_points=length(G_bay(:,1));
nr_loma_points=length(G_lomaprieta(:,1));

[x_channels, y_channels]=getCoordinates;

x_channels=x_channels(channel_indexes); y_channels=y_channels(channel_indexes);
x_channels=x_channels(compare_elements); y_channels=y_channels(compare_elements);



%% Compute predicted topography and quantiles FIX!!!

% load(filein,'x')

elev_pred = nan(length(x_channels)+nr_bay_points+nr_loma_points,round((length(shear))/100));
j = 0;

timerVal = tic;
for i = 1:100:length(shear)
% for i = 1:length(e11)
    j = j+1;
    mi = [shear(i); push(i); Kcon(:,i)];
    elev_pred(:,j) = elev_fun_ibem_n1_platemotion_Bay_Loma_multiK(mi,G_tect,G_bay,G_lomaprieta,geo_map,d,Ginv_elev,channel_indexes,compare_elements);
end
elapsedTime_elev_pred = toc(timerVal);

elev_pred=elev_pred(1:length(compare_elements),:);
%%
timerVal = tic;
elev_quant = quantile(elev_pred',[0.025 0.5 0.975]);
elapsedTime_elev_pred_quant = toc(timerVal);
%%

elev_ML=elev_fun_ibem_n1_platemotion_Bay_Loma_multiK(m_ML,G_tect,G_bay,G_lomaprieta,geo_map,d,Ginv_elev,channel_indexes,compare_elements);
elev_ML=elev_ML(1:length(compare_elements));

nr_model=length(m_ML)-nr_nounits;
nr_data=length(elev_ML);

chi_sq=(1/((nr_data-nr_model)*sig_elev^2))*...
    (sum((channel_elevations_to_compare-elev_ML).^2));

corr_var=sqrt((sig_elev^2)*(chi_sq));

resi_elev_ML=channel_elevations_to_compare-elev_ML;

%%
col=ones(nr_data,3)./255;
% col(geo_map(i,2):geo_map(i,3),3).*geocol{6}(i)]./255;

for i=1:Ngeo-nr_nounits
    col(geo_map(i,2):geo_map(i,3),:) =...
        [col(geo_map(i,2):geo_map(i,3),1).*geocol{4}(i) ...
        col(geo_map(i,2):geo_map(i,3),2).*geocol{5}(i) ...
        col(geo_map(i,2):geo_map(i,3),3).*geocol{6}(i)];
end

%% Save reduced chi_squared and corrected variance

filename_model_redchisq_corrvar=[weight_flag,'_weight_ML_redchisq_corrvariance.txt'];


fileID = fopen(filename_model_redchisq_corrvar,'w');
fprintf(fileID,'%7.4f %7.4f\n',[chi_sq, corr_var]');
fclose(fileID);

%%
figure

subplot(1,3,1), hist(channel_elevations_to_compare,100);
title([title_analysis, '\newline Hist Elevations Data'])
xlim([0 1000])
xlabel('Elevation [m]')
% daspect([0.5 1 1])
% font(16)

subplot(1,3,2), hist(elev_ML,100);

title([title_analysis, '\newline Hist Maximum Likelihood Elevations'])
xlim([0 1000])
xlabel('Elevation [m]')
% daspect([0.5 1 1])
% font(16)

subplot(1,3,3), hist(resi_elev_ML,100)
title([title_analysis, '\newline Hist Residuals Data-ML'])
xlim([-500 500])
xlabel('Elevation [m]')
% daspect([0.1 1 1])
% font(16)

set(gcf,'Position',[0, 0, 1400, 350])
% set(gcf,'position',[150   372   900   300])

%% Paper ready
figure
histogram(resi_elev_ML,100,'FaceColor','k','EdgeColor',[.8 .8 .8])
% histogram(resi_elev_ML,100,'EdgeAlpha',0.5)

title([title_analysis, '\newline Hist Residuals Data-ML'])
xlim([-500 500])
xlabel('Elevation [m]')
font(16)
set(gcf,'Position',[0, 0, 500, 700])
%%
figure

subplot(1,2,1)
dscatter(channel_elevations_to_compare,elev_ML)
% colorbar
colormap(jet)
axis equal
xlim([0 1100])
ylim([0 1100])
label('Channel elevations data','Predicted elevations Max Like')

subplot(1,2,2)
dscatter(log10(channel_elevations_to_compare),log10(elev_ML))
% colorbar
colormap(jet)
axis equal
xlim([1.6 3.2])
ylim([1.6 3.2])
label('Log Channel elevations data','Log Predicted elevations Max Like')

set(gcf,'position',[150   372   750   350])
supertitle(['Dscatter plot Predicted vs Data ', title_analysis])
font(16)

set(gcf,'Position',[0, 0, 1200, 500])

%%
figure

subplot(1,2,1)
scatter(x_channels./1e3,y_channels./1e3,10,channel_elevations_to_compare,'filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Channel Elevations Data')
font(16)

subplot(1,2,2)
scatter(x_channels./1e3,y_channels./1e3,10,elev_ML,'filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title('Channel Elevations Maximum Likelihood')
font(16)

supertitle([title_analysis,'\newline '])
font(16)

set(gcf,'Position',[0, 0, 1400, 500])

%%
figure

subplot(1,3,1)
scatter(x_channels./1e3,y_channels./1e3,10,elev_quant(1,:)','filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Channel Elev Predicted Lower 2.5 Quantile')
font(16)

subplot(1,3,2)
scatter(x_channels./1e3,y_channels./1e3,10,elev_ML,'filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title('Channel Elev Maximum Likelihood')
font(16)

subplot(1,3,3)
scatter(x_channels./1e3,y_channels./1e3,10,elev_quant(3,:)','filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title('Channel Elev Predicted Upper 97.5 Quantile')
font(16)

supertitle([title_analysis,'\newline '])
font(16)

set(gcf,'Position',[0, 0, 1400, 350])

%%
figure

subplot(1,3,1)
scatter(x_channels./1e3,y_channels./1e3,10,channel_elevations_to_compare,'filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Channel Elevations Data')
font(16)

subplot(1,3,2)
scatter(x_channels./1e3,y_channels./1e3,10,elev_ML,'filled')
caxis([0 1000]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title('Channel Elev Maximum Likelihood')
font(16)

subplot(1,3,3)
scatter(x_channels./1e3,y_channels./1e3,10,resi_elev_ML,'filled')
caxis([-300 300]), colorbar, colormap(jet)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([550,640])
ylim([4075,4150])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title('Residual')
font(16)

supertitle([title_analysis,'\newline '])
font(16)

set(gcf,'Position',[0, 0, 1400, 350])


%%

figure
scatter(x_channels./1e3,y_channels./1e3,8,col,'filled')
% caxis([0 1000]), colorbar, colormap(jet)
% hold on
% plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title('Channel Elevations Data [m]')
font(16)
box on
% hold off
set(gcf,'Position',[0, 0, 900, 750])
% set(gcf,'Position',[0, 0, 800, 600])


%% Paper ready
ax=figure;
scatter(x_channels./1e3,y_channels./1e3,8,channel_elevations_to_compare,'filled')
caxis([0 1000])
% colorbar, colormap(parula)
colorbar, cptcmap('GMT_drywet',ax,'flip',true)
% hold on
% plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title(['Channel Elevations Data [m]'])
font(16)
box on
% hold off
set(gcf,'Position',[0, 0, 900, 750])


ax=figure;
scatter(x_channels./1e3,y_channels./1e3,10,elev_ML,'filled')
caxis([0 1000])
% colorbar, colormap(jet)
colorbar, cptcmap('GMT_drywet',ax,'flip',true)
hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title([title_analysis,', Channel Elev Maximum Likelihood [m]'])
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])


%% Paper ready
ax=figure;
scatter(x_channels./1e3,y_channels./1e3,8,resi_elev_ML,'filled')
caxis([-300 300])
% colorbar, colormap(ax,redblue)
% colorbar, cptcmap('GMT_polar',ax,'flip',false)
CT=cbrewer('div', 'RdGy', 100); colorbar, colormap(ax,CT)

hold on
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([560,640])
ylim([4080,4150])
xlabel('UTM E [km]')
ylabel('UTM N [km]')
% title(['Channel Maximum Likelihood\newline ',title_analysis])
title([title_analysis,', Residual (Data-ML) [m]'])
font(16)
box on
hold off
set(gcf,'Position',[0, 0, 900, 750])


%%
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries')
addpath(genpath('/Users/felipearon/Documents/POSTDOC/Codes/BEM/Jacks_tribemx/tribem2018'))
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh')
addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/CEES_runs/')

% chi_file='obs_data_chi.txt';
grid_file='obs_data_grid.txt';

MCMC_MaxLike_shear_push=m_ML(1:2)*1e3;  %in mm/yr

cl='2.5';
% v=3;

% xyz_chi=load(chi_file);
xyz_chi=[x_channels,y_channels,channel_elevations_to_compare];
xyz_grid=load(grid_file);

len_chi=length(channel_elevations_to_compare);
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

faults=ReadPatches(faultnames);

ends = cumsum(faults.nEl); % Give ending indices of the 3 faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the 3 faults

inname_shear=['GreenFncs_cl',cl,'_shear.mat']; load(inname_shear);
inname_push=['GreenFncs_cl',cl,'_push.mat']; load(inname_push);

%%%%%FIX THIS!%%%%%
u_shear=reshape(obs_shear.u,3,[])';
u_shear_chi=u_shear(1:len_chi,3);
u_shear_grid=u_shear(len_chi+1:end,3);

u_push=reshape(obs_push.u,3,[])';
u_push_chi=u_push(1:len_chi,3);
u_push_grid=u_push(len_chi+1:end,3);
%%%%%%%%%%%%%%%%%%%%
%%
str_model=[slip_shear(:,1) slip_push(:,1)]*(MCMC_MaxLike_shear_push);
dip_model=[slip_shear(:,2) slip_push(:,2)]*(MCMC_MaxLike_shear_push);
nor_model=[slip_shear(:,3) slip_push(:,3)]*(MCMC_MaxLike_shear_push);

ux_model=[u_shear(:,1) u_push(:,1)]*(MCMC_MaxLike_shear_push);
ux_model_chi=ux_model(1:len_chi);
ux_model_grid=ux_model(end-len_grid+1:end);

uy_model=[u_shear(:,2) u_push(:,2)]*(MCMC_MaxLike_shear_push);
uy_model_chi=uy_model(1:len_chi);
uy_model_grid=uy_model(end-len_grid+1:end);

uz_model=[u_shear(:,3) u_push(:,3)]*(MCMC_MaxLike_shear_push);
uz_model_chi=uz_model(1:len_chi);
uz_model_grid=uz_model(end-len_grid+1:end);


%%
figure
% quiver(xyz_grid(:,1),xyz_grid(:,2),u_shear(len_chi+1:end,1),u_shear(len_chi+1:end,2),1,'-mo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
quiver(xyz_grid(:,1),xyz_grid(:,2),ux_model_grid,uy_model_grid,1,'m')
hold on
plot(xyz_grid(:,1),xyz_grid(:,2),'og','MarkerFaceColor','k','MarkerSize',1.5)
plot(traces(:,1),traces(:,2),'-k')
axis equal
xlim([380,830])
ylim([3880,4355])
hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title(['Grid xy-disp Maximum Likelihood\newline ',title_analysis])
set(gcf,'Position',[0, 0, 800, 600])
font(16)

%%
ax=figure;
scatter(xyz_chi(:,1)./1e3,xyz_chi(:,2)./1e3,8,uz_model_chi,'filled')
% hold on
% plot(traces(:,1),traces(:,2),'-k')
caxis([0 max(uz_model_chi)])
% colorbar, colormap(jet)
CT=cbrewer('seq', 'OrRd', 100); colorbar, colormap(ax,CT)
axis equal
xlim([560,640])
ylim([4080,4150])
% hold off
xlabel('UTM E [km]')
ylabel('UTM N [km]')
title(['Uplift Rate [mm/y] chi obs points Maximum Likelihood\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
font(16)
box on
% hold off
set(gcf,'Position',[0, 0, 900, 750])



%%
meshview(faults.c, faults.v(begs(1):ends(5),:));
% colorbar, colormap(redblue), caxis([min(str_MCMC(begs(3):ends(3))), min(str_MCMC(begs(3):ends(3)))*-1])
axis equal
% view(-35,0)
% shading flat
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title('Geometric Arrangment')
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
% set(gcf,'Position',[0, 0, 1200, 180])
%%
meshview(faults.c, faults.v(begs(3):ends(3),:), str_model(begs(3):ends(3)));
colorbar, colormap(redblue), caxis([min(str_model(begs(3):ends(3))), min(str_model(begs(3):ends(3)))*-1])
axis equal
view(-35,0)
shading flat
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Strike-slip [mm/yr] Max Like San Andreas (-)right-lat, (+)left-lat\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 180])

meshview(faults.c, faults.v(begs(3):ends(3),:), str_model(begs(3):ends(3)));
colorbar, colormap(redblue), caxis([min(str_model(begs(3):ends(3))), min(str_model(begs(3):ends(3)))*-1])
axis equal
view(-35,0)
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Strike-slip [mm/yr] Max Like San Andreas (-)right-lat, (+)left-lat\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 180])

%%

meshview(faults.c, faults.v(begs(1):ends(1),:), dip_model(begs(1):ends(1)));
colorbar, colormap(redblue), caxis([-1*max(dip_model(begs(2):ends(2))) max(dip_model(begs(2):ends(2)))]),
axis equal
shading flat
view(-7,15)
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Dip-slip [mm/yr] Max Like Berrocal (-)normal, (+)thrust\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 300])

meshview(faults.c, faults.v(begs(1):ends(1),:), dip_model(begs(1):ends(1)));
colorbar, colormap(redblue), caxis([-1*max(dip_model(begs(2):ends(2))) max(dip_model(begs(2):ends(2)))]),
axis equal
view(-7,15)
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Dip-slip [mm/yr] Max Like Berrocal (-)normal, (+)thrust\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 300])

%%
meshview(faults.c, faults.v(begs(2):ends(2),:), dip_model(begs(2):ends(2)));
colorbar, colormap(redblue), caxis([-1*max(dip_model(begs(2):ends(2))) max(dip_model(begs(2):ends(2)))]),
axis equal
shading flat
view(-7,15)
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Dip-slip Max Like [mm/yr] Shannon-Mt Vista (-)normal, (+)thrust\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 300])

meshview(faults.c, faults.v(begs(2):ends(2),:), dip_model(begs(2):ends(2)));
colorbar, colormap(redblue), caxis([-1*max(dip_model(begs(2):ends(2))) max(dip_model(begs(2):ends(2)))]),
axis equal
% shading flat
view(-7,15)
xlabel('UTM E [km]')
ylabel('UTM N [km]')
zlabel('Depth [km]')
title(['Dip-slip Max Like [mm/yr] Shannon-Mt Vista (-)normal, (+)thrust\newline ',title_analysis])
% title('Uz (uplift) chi obs points MCMC MaxLike finest model (18,621 elements)')
% set(gcf,'Position',[0, 0, 800, 600])
font(16)
set(gcf,'Position',[0, 0, 1200, 300])

% vline(q(1),'--r',2);
% vline(q(2),'--r',2);
% vline(e11stats.ML*1e9,'-w',6);
% vline(e11stats.ML*1e9,'-r',2);
% vline(e22stats.ML*1e9,'-w',6);
% vline(e22stats.ML*1e9,'-r',2);

%%%%finish edits%%%%











%% Plot K values
keys = textread('geo_keys.txt','%s');
fid = fopen('geokeys.txt');
geocol = textscan(fid,'%s%d%d%f%f%f%s%d');
fclose(fid);

% Kquant = Kstats.Kquant;
Kxx = Kstats.Kxx;
Kf = Kstats.Kf;
% Kml = Kstats.Kml;

plotorder = double(geocol{8});

K_ML = m_ML(4:end);
    
figure
for k = 1:Ngeo
    i = find(ismember(plotorder,k));
    subplot(5,4,k); hold on
    col = [geocol{4}(i) geocol{5}(i) geocol{6}(i)]./255;
    area(Kxx(i,:),Kf(i,:),'facecolor',col,'linewidth',2)
    %     title(['log K_{' num2str(i) '}'])
    title(keys(i))
    font(16)
    
    %     vline(Kml(i),'-r',2);
    %     vline(Kquant(i,1),'--r',2);
    %     vline(Kquant(i,2),'--r',2);
    
    vline(K_ML(i),'--r',2);
    
    
    if k == 1
        label('log K', 'Prob. density',16)
    end
    axis tight
    xlim([-6 -1])

    
    box on
    
end


%% Plot K individually 

for k = 1:Ngeo
    figure
    i = find(ismember(plotorder,k));
    hold on
    col = [geocol{4}(i) geocol{5}(i) geocol{6}(i)]./255;
    area(Kxx(i,:),Kf(i,:),'facecolor',col,'linewidth',2)
    %     title(['log K_{' num2str(i) '}'])
    title(keys(i))
    font(16)
    
    %     vline(Kml(i),'-r',2);
    %     vline(Kquant(i,1),'--r',2);
    %     vline(Kquant(i,2),'--r',2);
    
    vline(K_ML(i),'--r',2);
    
    
    label('log K', 'Prob. density',16)
    axis tight
    xlim([-6 -1])

    
    box on
    
end



%% Compute predicted Ksn and quantiles

% load(filein,'x')



Ksn_pred = nan(length(x),round((length(e11))/100));
j = 0;
for i = 1:100:length(e11)
% for i = 1:length(e11)
    j = j+1;
    mi = [e11(i); e22(i); e12(i); Kcon(:,i)]; 
    Ksn_pred(:,j) = Ksn_fun_ibem(mi,G,geo);
end

Ksn_quant = quantile(Ksn_pred',[0.025 0.5 0.975]);


Ksn_ML = Ksn_fun_ibem(m_ML,G,geo);

%% Plot predicted Ksn, compare with data

figure

% Lower 95% confidence bound
% subplot(231)
scatter(x,y,20,Ksn_quant(1,:),'filled')
colorbar
caxis([25 300])
axis equal
title('Lower 95% confidence bound');
colormap(jet)

% Max lik
% subplot(232)
figure
scatter(x,y,20,Ksn_ML,'filled')
colorbar
caxis([25 300])
axis equal
title('Maximum likelihood');
colormap(jet)

% Upper 95% confidence bound
% subplot(233)
figure
scatter(x,y,20,Ksn_quant(3,:),'filled')
colorbar
caxis([25 300])
axis equal
title('Upper 95% confidence bound');
colormap(jet)

% subplot(235)
figure
scatter(x,y,20,ksn,'filled')
colorbar
caxis([25 300])
axis equal
title('Data');
colormap(jet)

% residual
res = Ksn_ML - ksn;
figure
scatter(x,y,20,res,'filled')
colorbar
caxis([-100 100])
axis equal
title('Residual');
colormap(jet)


%% Predicted K
% 
% K_pred = nan(length(x),length(Xkeep));
% for i = 1:length(Xkeep)
%     K_pred(:,i) = U./Ksn_pred(:,i);
% end
% 
% %
% % PK1 = U'./P(1,:);
% % % PK2 = U'./mean(Ksn_pred');
% % PK3 = U'./P(3,:);
% 
% PK = quantile(K_pred', [0.025 0.5 0.975]);
% 
% figure
% plot(x_obs,K,'r.','linewidth',2,'markersize',20)
% hold on
% % plot(x_obs,Ksn_all,'--')
% plot(x_obs,PK(1,:),'bo','linewidth',2)
% % plot(x_obs,PK2,'bo')
% plot(x_obs,PK(3,:),'bo','linewidth',2)
% label('x','K')
% grid
% legend('True','95% C.I.')


%% Plot uplift

% imax = find(logL_keep==max(logL_keep),1);

% upl = G*[sh_stats.ML; pu_stats.ML];
% upl = G*[e11(iML); e22(iML); e12(iML)]*1e6;  % in mm/yr
upl = G*[e11(iML); e22(iML); e12(iML)]*1e3;  % in mm/yr

figure
scatter(x,y,[],upl,'filled')
colorbar
title('uplift (mm/yr)')
axis equal
ax = axis;

% caxis([0 3])

load FaultTraces
hold on
plot(BerrocalTrace(:,1)*1e3,BerrocalTrace(:,2)*1e3,'k')
plot(ShannonTrace(:,1)*1e3,ShannonTrace(:,2)*1e3,'k')
plot(SanAndreasTrace(:,1)*1e3,SanAndreasTrace(:,2)*1e3,'k')

axis(1e6*[0.5511    0.6480    4.0665    4.1559])

colormap(jet)

% axis(ax);

%% Plot map of estimated Ks


geo_codes = unique(geo);
Ngeo = length(geo_codes);

% % IND = cell(Ngeo,1);
% figure
% hold on
% for i = 1:length(geo_codes)
%     IND = find(geo == geo_codes(i));   
%     scatter(x(IND),y(IND),[],Kstats.Kquant(i,1)*ones(length(IND),1),'filled')
% end
% axis equal
% colorbar
% caxis([-8.4 -8])
% title('log(K) Lower 95% bound')

figure
hold on
for i = 1:length(geo_codes)
    IND = find(geo == geo_codes(i));   
    scatter(x(IND),y(IND),[],Kcon(i,iML)*ones(length(IND),1),'filled')
end
axis equal
colorbar
caxis([-6 -4.8])
title('log(K) Maximum likelihood')
% 
% figure
% hold on
% for i = 1:length(geo_codes)
%     IND = find(geo == geo_codes(i));   
%     scatter(x(IND),y(IND),[],Kstats.Kquant(i,2)*ones(length(IND),1),'filled')
% end
% axis equal
% colorbar
% caxis([-7.7 -7.2])
% title('log(K) Upper 95% bound')


load FaultTraces
hold on
plot(BerrocalTrace(:,1)*1e3,BerrocalTrace(:,2)*1e3,'k')
plot(ShannonTrace(:,1)*1e3,ShannonTrace(:,2)*1e3,'k')
plot(SanAndreasTrace(:,1)*1e3,SanAndreasTrace(:,2)*1e3,'k')
axis(1e6*[0.5511    0.6480    4.0665    4.1559])


%%

% Kml = 10.^Kcon(:,iML)
% ksn_for = 1e-3*upl./(;


%% K stats for each observation point

% K_ML = nan(length(x),1);
% K_U95 = K_ML;
% K_L95 = K_ML;
% 
% 
% for i = 1:length(geo_codes)
%     IND = find(geo == geo_codes(i)); 
%     K_ML(IND) = Kcon(i,iML)*ones(length(IND),1);
% %     K_U95(IND) = Kstats.Kquant(i,2)*ones(length(IND),1);
% %     K_L95(IND) = Kstats.Kquant(i,1)*ones(length(IND),1);
% 
% %     scatter(x(IND),y(IND),[],Kcon(i,iML)*ones(length(IND),1),'filled')
% end
% 
% KSTATS.K_ML = K_ML;
% KSTATS.K_U95 = K_U95;
% KSTATS.K_L95 = K_L95;




%% Plot correlations

figure;

subplot(221)
hold on
dscatter(-e11(1:100:end)'*1e9,e12(1:100:end)'*1e9,'plottype','contour')
plot(-e11stats.ML*1e9,e12stats.ML*1e9,'pk','markersize',20,'markerfacecolor','c')
colormap(flipud(hot))
label('\epsilon_{11} (10^{-9})','\epsilon_{12} (10^{-9})',16)
axis equal
font(16)


subplot(222)
hold on
dscatter(-e22(1:100:end)'*1e9,e12(1:100:end)'*1e9,'plottype','contour')
plot(-e22stats.ML*1e9,e12stats.ML*1e9,'pk','markersize',20,'markerfacecolor','c')
colormap(flipud(hot))
label('\epsilon_{22} (10^{-9})','\epsilon_{12} (10^{-9})',16)
axis equal
font(16)

subplot(224)
hold on
dscatter(-e22(1:100:end)'*1e9,-e11(1:100:end)'*1e9,'plottype','contour')
plot(-e22stats.ML*1e9,-e11stats.ML*1e9,'pk','markersize',20,'markerfacecolor','c')
colormap(flipud(hot))
label('\epsilon_{22} (10^{-9})','\epsilon_{11} (10^{-9})',16)
axis equal
font(16)
grid on