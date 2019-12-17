scenario = 'strong';

% Setting the correct directories
p=pathdef; path(p)

load(strcat('results/multiK_', scenario, '_enforcement_v_with_cov_linearized'));

vs_mean = mean(v(:,1))*1000;
vn_mean = mean(v(:,2))*1000;
vs_std = std(v(:,1))*1000;
vn_std = std(v(:,2))*1000;

fprintf('Mean vs (mm/yr): %6.2f', vs_mean);
fprintf('Mean vn (mm/yr): %6.2f', vn_mean);
fprintf('Std vs (mm/yr): %6.2f', vs_std);
fprintf('Std vn (mm/yr): %6.2f', vn_std);

x_s = vs_mean-3*vs_std:0.01:vs_mean+3*vs_std;
x_n = vn_mean-3*vn_std:0.01:vn_mean+3*vn_std;

[X_S, X_N] = meshgrid(x_s, x_n);

pdf_s = pdf(v(:,1)*1000, x_s);
pdf_n = pdf(v(:,2)*1000, x_n);

pdf_sn = pdf([v(:,1:2)*1000], [X_S, X_N]);

figure(1)
plot(x_s, pdf_s, 'k-');
xlabel('Shear Velocity (mm/yr)');
ylabel('Density (yr/mm)');

figure(2)
plot(x_n, pdf_n, 'k-');
xlabel('Normal Velocity (mm/yr)');
ylabel('Density (yr/mm)');

figure(3)
imagesc(X_S, X_N, pdf_sn);
xlabel('Shear Velocity (mm/yr)');
ylabel('Normal Velocity (mm/yr)');
colorbar;
