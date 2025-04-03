close all

wt_hb = [484 57 98; 494 102 90.15; 493 101 22.19; 449 100 90.5; 493 104 53; 490 104 37.66 ];
wt_hb = [wt_hb(:,1)-333 wt_hb(:,2) wt_hb(:,3)];

wt_extract_rbd = [453 495 483 486 488 451 491 ];
wt_extract_nb = [32 58 59];
% Add hydrogen bond percentages
for i = 1:size(wt_hb,1)
    wt_dist_cov(wt_hb(i,1),wt_hb(i,2))=wt_hb(i,3);
end

% Exclusions
for i = 1:size(wt_extract_rbd,2)
    wt_dist_cov(wt_extract_rbd(i)-333,:)=0;
end

% Exclusions
for i = 1:size(wt_extract_nb,2)
    wt_dist_cov(:,wt_extract_nb(i))=0;
end

wt_dist_cov(wt_dist(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
figure(2)
subplot(2,2,1)

h1 = imagesc(rbd_int, nb_ind, wt_dist_cov(rbd_int-333,nb_ind)');
set(gca,'YDir','normal')
grid off
hold on
colorbar

% % Define the y-values for the horizontal lines
% y_values = [26 32 52 56 99 116];
% 
% % Loop to plot the horizontal lines
% for y = y_values
%     plot([445 506], [y y], "-k")
% end

xlim([445 500])
% ylim([20 120])
xlabel("RBD Residue")
ylabel("Nanobody Residues")
title("WT")
clim([-100 100])
cb=colorbar;
cmap=colormap;
cmap(1:127,:) = [linspace(0,1,127); linspace(0,1,127); linspace(1,1,127)]';
cmap(128:129,:)=[1 1 1; 1 1 1];
cmap(130:256,:)=[linspace(1,1,127); linspace(1,0,127); linspace(1,0,127)]';
ax= gca;
ax.Colormap=cmap;


subplot(2,2,2)

h2 = imagesc(rbd_int, nb_ind, alpha_dist_cov(rbd_int-333,nb_ind)');
set(gca,'YDir','normal')
grid off
hold on
colorbar

% % Define the y-values for the horizontal lines
% y_values = [26 32 52 56 99 116];
% 
% % Loop to plot the horizontal lines
% for y = y_values
%     plot([445 506], [y y], "-k")
% end

xlim([445 500])
% ylim([20 120])
xlabel("RBD Residue")
ylabel("Nanobody Residues")
title("Alpha")
clim([-100 100])
cb=colorbar;
cmap=colormap;
cmap(1:127,:) = [linspace(0,1,127); linspace(0,1,127); linspace(1,1,127)]';
cmap(128:129,:)=[1 1 1; 1 1 1];
cmap(130:256,:)=[linspace(1,1,127); linspace(1,0,127); linspace(1,0,127)]';
ax= gca;
ax.Colormap=cmap;


subplot(2,2,3)

h3 = imagesc(rbd_ind, nb_ind, beta_dist_cov(rbd_ind-333,nb_ind)');
set(gca,'YDir','normal')
grid off
hold on
colorbar

% % Define the y-values for the horizontal lines
% y_values = [26 32 52 56 99 116];
% 
% % Loop to plot the horizontal lines
% for y = y_values
%     plot([445 506], [y y], "-k")
% end

xlim([445 500])
% ylim([20 120])
xlabel("RBD Residue")
ylabel("Nanobody Residues")
title("Beta")
clim([-100 100])
cb=colorbar;
cmap=colormap;
cmap(1:127,:) = [linspace(0,1,127); linspace(0,1,127); linspace(1,1,127)]';
cmap(128:129,:)=[1 1 1; 1 1 1];
cmap(130:256,:)=[linspace(1,1,127); linspace(1,0,127); linspace(1,0,127)]';
ax= gca;
ax.Colormap=cmap;

subplot(2,2,4)

h4 = imagesc(rbd_int, nb_ind, omicron_dist_cov(rbd_int-333,nb_ind)');
set(gca,'YDir','normal')
grid off
hold on
colorbar

% % Define the y-values for the horizontal lines
% y_values = [26 32 52 56 99 116];
% 
% % Loop to plot the horizontal lines
% for y = y_values
%     plot([445 506], [y y], "-k")
% end

xlim([445 500])
% ylim([20 120])
xlabel("RBD Residue")
ylabel("Nanobody Residues")
title("Omicron")
clim([-100 100])
cb=colorbar;
cmap=colormap;
cmap(1:127,:) = [linspace(0,1,127); linspace(0,1,127); linspace(1,1,127)]';
cmap(128:129,:)=[1 1 1; 1 1 1];
cmap(130:256,:)=[linspace(1,1,127); linspace(1,0,127); linspace(1,0,127)]';
ax= gca;
ax.Colormap=cmap;


figure(3)
h1 = imagesc(rbd_int, nb_ind, wt_dist_cov(rbd_int-333,nb_ind)');
set(gca,'YDir','normal')
grid off
hold on
colorbar

% % Define the y-values for the horizontal lines
% y_values = [26 32 52 56 99 116];
% 
% % Loop to plot the horizontal lines
% for y = y_values
%     plot([445 506], [y y], "-k")
% end

% xlim([445 500])
% % ylim([20 120])
xlabel("RBD Residue")
ylabel("Nanobody Residues")
title("WT")
clim([-100 100])
cb=colorbar;
cmap=colormap;
cmap(1:127,:) = [linspace(0,1,127); linspace(0,1,127); linspace(1,1,127)]';
cmap(128:129,:)=[1 1 1; 1 1 1];
cmap(130:256,:)=[linspace(1,1,127); linspace(1,0,127); linspace(1,0,127)]';
ax= gca;
ax.Colormap=cmap;