clear all;
FOV_X = 20; N = 100;
in_plane_coord = linspace(-FOV_X/2,FOV_X/2,N);  % coordinates for output field map
[YY,ZZ] = meshgrid(in_plane_coord,in_plane_coord);


figure;
xslicevec = -100:10:100;
for ii = 1:numel(xslicevec)
xslice = xslicevec(ii);
load(['FM',num2str(xslice),'mm_4_18_17.mat']);
F = TriScatteredInterp(pos(:,1),pos(:,2), fm_freq.');
vq(:,:,ii) = F(YY,ZZ)*42.576e3;
subplot(5,5,ii)
mesh(YY,ZZ,vq(:,:,ii));
% hold on; plot3(pos(:,1),pos(:,2), fm_freq.','o'); hold off; pause(0.1);
title(['x = ',num2str(xslice)]); colormap jet; axis square;
%  axis([0 1000 -0.001 110 0 110])    
 view(0,90); caxis([-40e3,40e3]);
end

fm_freq = fm_freq*42.576e3;

I=[];
[I,Y]  = find(pos(:,2) == 0);
[Y,I2]  = max(pos(I,1));
I3 = I(I2);
I=[];
[I,Y]  = find(pos(:,2) == 0);
[Y,I2]  = min(pos(I,1));
I4 = I(I2);
fm_freq(I4) - fm_freq(I3)


I=[];
[I,Y]  = find(pos(:,1) == 0);
[Y,I2]  = max(pos(I,2));
I3 = I(I2);
I=[];
[I,Y]  = find(pos(:,1) == 0);
[Y,I2]  = min(pos(I,2));
I4 = I(I2);
fm_freq(I4) - fm_freq(I3)






pos_all = [];  B_all = [];
for ii = -100:10:100

load(['FM',num2str(ii),'mm_4_18_17.mat']);
temp = ones(size(pos(:,1))).*ii/10;
pos = [temp, pos];
pos_all = [pos_all; pos];
B_all = [B_all, fm_freq];
end




orders_to_calculate = [0:5];
unit_sph_harm_datapoints = calc_spherical_harmonics_arb_points(orders_to_calculate,pos_all);


%% Solve for spherical harmonic coefficients based on data
coeff = unit_sph_harm_datapoints\B_all.';

%% generate field maps for output =========================================
% generate coordinate vector for each slice
FOV_X = 20; N = 41;
in_plane_coord = linspace(-FOV_X/2,FOV_X/2,N);  % coordinates for output field map

% [YY,ZZ] = ndgrid(in_plane_coord,in_plane_coord);
% YY_vec = reshape(YY,numel(YY),1);  ZZ_vec = reshape(ZZ,numel(ZZ),1);


% xvec = [-10:0.5:10];  N=numel(xvec);
% pos_cube = ptgrid_cube(xvec).';

[X,Y,Z] = meshgrid(in_plane_coord,in_plane_coord,in_plane_coord);
pos_all = [X(:),Y(:),Z(:)].';

% load unit_sph_harm_8order_5mmres;  
% load unit_sph_harm_5order_5mmres_cube;
% 
% % slice_centers = 0; ss = 1;
%    
%     % this code takes a while for large grid sizes 'N' above ~64... 
%      temp = calc_spherical_harmonics_arb_points(orders_to_calculate,[xslice*ones(numel(YY_vec),1) YY_vec ZZ_vec]);
% unit_sph_harm = calc_spherical_harmonics_arb_points(orders_to_calculate,pos_cube.');
%   load('unit_sph_harm_5order_5mmres_cmunits');
 
unit_sph_harm = calc_spherical_harmonics_arb_points(orders_to_calculate,pos_all.');

total_field = unit_sph_harm*coeff*42.576e3;
% fieldmap = reshape(total_field,N,N,N)*42.58e3;

temp = reshape(total_field,N,N,N);
% clim = [min(temp(:)),max(temp(:))];
clim = [-40e3,40e3];
figure;
subplot(1,3,1)
imagesc(linspace(-10,10,200), linspace(-10,10,200), rot90(squeeze(temp(:,:,ceil(end/2))))); axis square; colormap jet; caxis(clim);
xlabel('X (cm)'); ylabel('Y (cm)');
subplot(1,3,2)
imagesc(linspace(-10,10,200), linspace(-10,10,200), rot90(squeeze(temp(ceil(end/2),:,:)))); axis square; colormap jet; caxis(clim);
xlabel('Y (cm)'); ylabel('Z (cm)');
subplot(1,3,3)
imagesc(linspace(-10,10,200), linspace(-10,10,200), rot90(squeeze(temp(:,ceil(end/2),:)))); axis square; colormap jet; caxis(clim);
xlabel('X (cm)'); ylabel('Z (cm)');
% fieldmap = reshape(total_field,N,N)*42.58e6;
%     
% fave = mean(mean(fieldmap));
% 
% figure;
% imagesc(linspace(-10,10,200), linspace(-10,10,200), rot90(squeeze(fieldmap(ceil(end/2),:,:)))); axis square; %title(['xslice =  ',num2str(xslice),'cm']); 
% colormap jet;
% % 
% % 
% % caxis([fave-30e3, fave+30e3]);
% %  
% %  hold on;
% %  
% % 
% %  circle(0,0,sqrt(10^2-abs(xslice)^2));
% % 
