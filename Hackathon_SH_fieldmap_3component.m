
clear all;
% data = readmatrix('cube_2p5cm_measurement.csv');
load('fieldmap_raw_231017.mat');
r =50;

index = data(:,1);
dx = data(:,2);
dy = data(:,3);
dz = data(:,4);
Bx = data(:,5)/1e4;
By = data(:,6)/1e4;
Bz = data(:,7)/1e4;

traj(1,:) = [0,0,0];

for nn = 2:numel(index)
    traj(nn,1) = traj(nn-1,1)+dx(nn);
    traj(nn,2) = traj(nn-1,2)+dy(nn);
    traj(nn,3) = traj(nn-1,3)+dz(nn);
end

pos(:,1) = traj(:,3);
pos(:,2) = traj(:,2);
pos(:,3) = traj(:,1);
field_matrix(:,1) = Bx;
field_matrix(:,2) = By;
field_matrix(:,3) = Bz;
field_matrix(:,4) = sqrt(Bx.^2+By.^2+Bz.^2);


field_map_order = 7;
% data_use = Bx;%field_matrix(:,4);
orders_to_calculate = [0:field_map_order];
unit_sph_harm_FM = calc_spherical_harmonics_arb_points(orders_to_calculate,pos);


% data_use = Bx;%field_matrix(:,4);
sph_coeffs_Bx = unit_sph_harm_FM\Bx;
sph_coeffs_By = unit_sph_harm_FM\By;
sph_coeffs_Bz = unit_sph_harm_FM\Bz;


res = 2.5;
xyi = -50:res:50;
N = numel(xyi);
[Yq,Xq,Zq] = meshgrid(xyi);

res = 2.5;
xyi = -50:res:50;
N = numel(xyi);
pos_fit = ptgrid_cube(xyi);

orders_to_calc = [0:field_map_order];
%     unit_sph_harm = calc_spherical_harmonics_arb_points(orders_to_calc,pos_fit);
%     save(['unit_sph_harm_',num2str(field_map_order),'order_',num2str(res),'mmres_cube_2023'],'unit_sph_harm');
load('unit_sph_harm_7order_2.5mmres_cube_2023.mat');

total_field_Bx = unit_sph_harm*sph_coeffs_Bx;
total_field_By = unit_sph_harm*sph_coeffs_By;
total_field_Bz = unit_sph_harm*sph_coeffs_Bz;

field_mat_Bx = reshape(total_field_Bx,N,N,N);
field_mat_By = reshape(total_field_By,N,N,N);
field_mat_Bz = reshape(total_field_Bz,N,N,N);

% Rq = sqrt(Xq.^2+Yq.^2+Zq.^2);
Rq = sqrt(pos_fit(:,1).^2+pos_fit(:,2).^2+pos_fit(:,3).^2);
[Iq] = find(Rq <= r);

field_sphere_Bx = total_field_Bx(Iq);
field_sphere_By = total_field_By(Iq);
field_sphere_Bz = total_field_Bz(Iq);

field_sphere_Bx2 = field_mat_Bx(Iq);
field_sphere_By2 = field_mat_By(Iq);
field_sphere_Bz2 = field_mat_Bz(Iq);


%% homogeneity reporting
Brangeq = range(field_sphere_Bx);
Bstdq = std(field_sphere_Bx);
Bmeanq = mean(field_sphere_Bx);
Bmean = mean(field_sphere_Bx);

display(['10cm DSV interp: B0 mean = ',num2str(Bmeanq,'%.2f'),'mT, range = ',num2str(Brangeq,'%.2f'),'mT, std = ', num2str(Bstdq,'%.2f'),'mT']);
display(['Brange/Bmean ppm = ',num2str(1e6*Brangeq/Bmeanq,'%.2f')]);



%% plotting
temp = field_matrix(:,4);
meanB = mean(temp(:));
C1 = (meanB - 0.1);
C2 = (meanB + 0.1);

r = 50;
x=0;
y=0;
th = 0:pi/100:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

figure;

subplot(1,3,1);
imagesc([-50,50],[-50,50],rot90(squeeze(field_mat_Bx(:,ceil(end/2),:)))); axis square;
% caxis([Bmean-cxx,Bmean+cxx]);
caxis([C1,C2]);

colormap jet;
xlabel('x (m)');
ylabel('z (m)');
colorbar;
title('B0mag (mT) at y = 0');
hold on; plot(xunit, yunit, '-k');
axis square;


subplot(1,3,2);
imagesc([-50,50],[-50,50],rot90(squeeze(field_mat_Bx(ceil(end/2),:,:)))); axis square;
caxis([C1,C2]);
colormap jet;
xlabel('y (m)');
ylabel('z (m)');
colorbar;
title('B0mag (mT) at x = 0');
hold on; plot(xunit, yunit, '-k');
axis square;

subplot(1,3,3);
imagesc([-50,50],[-50,50],rot90(squeeze(field_mat_Bx(:,:,ceil(end/2))))); axis square;
caxis([C1,C2]);
colormap jet;
xlabel('x (m)');
ylabel('y (m)');
colorbar;
title('B0mag (mT) at z = 0');
hold on; plot(xunit, yunit, '-k');
axis square;
