scenario = 'strong';

% Setting the correct directories
p=pathdef; path(p)

load(strcat('results/multiK_', scenario, '_enforcement_v_with_cov_linearized'));

vs_mean = mean(v(:,1))*1000;
vn_mean = mean(v(:,2))*1000;
vs_std = std(v(:,1))*1000;
vn_std = std(v(:,2))*1000;

fprintf('Mean vs (mm/yr): %6.2f\n', vs_mean);
fprintf('Mean vn (mm/yr): %6.2f\n', vn_mean);
fprintf('Std vs (mm/yr): %6.2f\n', vs_std);
fprintf('Std vn (mm/yr): %6.2f\n', vn_std);

x_s = vs_mean-3*vs_std:(6*vs_std)/200:vs_mean+3*vs_std;
x_n = vn_mean-3*vn_std:(6*vs_std)/200:vn_mean+3*vn_std;

[X_S, X_N] = meshgrid(x_s, x_n);

pdf_s = ksdensity(v(:,1)*1000, x_s);
pdf_n = ksdensity(v(:,2)*1000, x_n);

pdf_sn = reshape(ksdensity([v(:,1:2)*1000], [X_S(:) X_N(:)]),length(x_n),length(x_s));

figure(1)
plot(x_s, pdf_s, 'k-');
xlabel('Shear Velocity (mm/yr)');
ylabel('Density (yr/mm)');

figure(2)
plot(x_n, pdf_n, 'k-');
xlabel('Normal Velocity (mm/yr)');
ylabel('Density (yr/mm)');

figure(3)
imagesc(x_n,x_s,pdf_sn')
set(gca,'ydir','normal');
axis equal;axis image;colormap jet
colorbar;
ylabel('Shear Velocity (mm/yr)');
xlabel('Normal Velocity (mm/yr)');
