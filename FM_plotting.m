clear all;
load('fieldmap_fit_231017.mat');

field_mat_mag = sqrt(field_mat_Bx.^2+field_mat_By.^2+field_mat_Bz.^2);

%% B0 is in the X direction of the magnet
field_mat_use = field_mat_Bx;  %%choose which component to plot




%% plotting
meanB = mean(field_mat_use(:));
C1 = (meanB - 0.1);
C2 = (meanB + 0.1);

r = 50;  % FOV radius in mm
x=0;
y=0;
th = 0:pi/100:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

fitmaxpos = max(pos_fit(:));

figure;
subplot(1,3,1);
imagesc([-fitmaxpos,fitmaxpos],[-fitmaxpos,fitmaxpos],rot90(squeeze(field_mat_use(:,ceil(end/2),:)))); axis square;
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
imagesc([-fitmaxpos,fitmaxpos],[-fitmaxpos,fitmaxpos],rot90(squeeze(field_mat_use(ceil(end/2),:,:)))); axis square;
caxis([C1,C2]);
colormap jet;
xlabel('y (m)');
ylabel('z (m)');
colorbar;
title('B0mag (mT) at x = 0');
hold on; plot(xunit, yunit, '-k');
axis square;

subplot(1,3,3);
imagesc([-fitmaxpos,fitmaxpos],[-fitmaxpos,fitmaxpos],rot90(squeeze(field_mat_use(:,:,ceil(end/2))))); axis square;
caxis([C1,C2]);
colormap jet;
xlabel('x (m)');
ylabel('y (m)');
colorbar;
title('B0mag (mT) at z = 0');
hold on; plot(xunit, yunit, '-k');
axis square;
