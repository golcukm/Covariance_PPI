% Title: DCD Data Reader and Covariance Matrix Constructor

% Authors: Mert Golcuk, PhD (Postdoctoral Researcher), Mert Gur, PhD (PI)

% Lab: Gur Lab, University of Pittsburgh School of Medicine
% Date: 2025

% Description:
% This script visualizes covariance matrices derived from molecular dynamics trajectory data.
% The resulting visualizations are useful for understanding structural dynamics and interactions in 
% computational biophysics studies.

% Requirements:


% Usage:
% - Specify the input DCD file.
% - Run the script to compute and output the covariance matrix.

% License:
% Available for academic and research purposes under MIT License.

close all

% Set variables for residue indices
rbd_ind = 334:528; % RBD residue indices
nb_ind = 1:128; % Nanobody residue indices

% Set variables for interaction indices
rbd_int = 445:506; % RBD interaction indices
nb_int = [26:32 52:56 99:116]; % Nanobody interaction indices

% Create a figure
figure(1)

% Plot for WT (Wild Type)
subplot(2,2,1)
h1 = imagesc(rbd_int, nb_ind, cov_wt_d(rbd_int-333,nb_ind)'); % Display covariance matrix for WT
set(gca,'YDir','normal') % Set Y-axis direction to normal
grid off % Disable grid
hold on % Hold the current plot
colorbar % Add colorbar

xlim([445 500]) % Set x-axis limits
xlabel("RBD Residue") % Label x-axis
ylabel("Nanobody Residues") % Label y-axis
title("WT") % Set title
clim([-5 1]) % Set color limits

% Define colormap for WT
cmap = colormap;
kp = 256/6; % Define keypoint for colormap
cmap(1:kp*5,:) = [linspace(0, 1, floor(kp*5)); linspace(0, 1, floor(kp*5)); linspace(1, 1, floor(kp*5))]';
cmap(floor(kp*5):floor(kp*5)+1,:) = [1 1 1; 1 1 1]; % White color for mid-range
cmap(floor(kp*5)+2:256,:) = [linspace(1, 1, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1))]';
ax = gca;
ax.Colormap = cmap; % Apply colormap

% Plot for Alpha variant
subplot(2,2,2)
h2 = imagesc(rbd_int, nb_ind, cov_alpha_d(rbd_int-333,nb_ind)'); % Display covariance matrix for Alpha
set(gca,'YDir','normal') % Set Y-axis direction to normal
grid off % Disable grid
hold on % Hold the current plot
colorbar % Add colorbar

xlim([445 500]) % Set x-axis limits
xlabel("RBD Residue") % Label x-axis
ylabel("Nanobody Residues") % Label y-axis
title("Alpha") % Set title
clim([-5 1]) % Set color limits

% Define colormap for Alpha
kp = 256/6; % Define keypoint for colormap
cmap(1:kp*5,:) = [linspace(0, 1, floor(kp*5)); linspace(0, 1, floor(kp*5)); linspace(1, 1, floor(kp*5))]';
cmap(floor(kp*5):floor(kp*5)+1,:) = [1 1 1; 1 1 1]; % White color for mid-range
cmap(floor(kp*5)+2:256,:) = [linspace(1, 1, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1))]';
ax = gca;
ax.Colormap = cmap; % Apply colormap

% Plot for Beta variant
subplot(2,2,3)
h3 = imagesc(rbd_int, nb_ind, cov_beta_d(rbd_int-333,nb_ind)'); % Display covariance matrix for Beta
set(gca,'YDir','normal') % Set Y-axis direction to normal
grid off % Disable grid
hold on % Hold the current plot
colorbar % Add colorbar

xlim([445 500]) % Set x-axis limits
xlabel("RBD Residue") % Label x-axis
ylabel("Nanobody Residues") % Label y-axis
title("Beta") % Set title
clim([-5 1]) % Set color limits

% Define colormap for Beta
kp = 256/6; % Define keypoint for colormap
cmap(1:kp*5,:) = [linspace(0, 1, floor(kp*5)); linspace(0, 1, floor(kp*5)); linspace(1, 1, floor(kp*5))]';
cmap(floor(kp*5):floor(kp*5)+1,:) = [1 1 1; 1 1 1]; % White color for mid-range
cmap(floor(kp*5)+2:256,:) = [linspace(1, 1, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1))]';
ax = gca;
ax.Colormap = cmap; % Apply colormap

% Plot for Omicron variant
subplot(2,2,4)
h4 = imagesc(rbd_int, nb_ind, cov_omicron_d(rbd_int-333,nb_ind)'); % Display covariance matrix for Omicron
set(gca,'YDir','normal') % Set Y-axis direction to normal
grid off % Disable grid
hold on % Hold the current plot
colorbar % Add colorbar

xlim([445 500]) % Set x-axis limits
xlabel("RBD Residue") % Label x-axis
ylabel("Nanobody Residues") % Label y-axis
title("Omicron") % Set title
clim([-5 1]) % Set color limits

% Define colormap for Omicron
kp = 256/6; % Define keypoint for colormap
cmap(1:kp*5,:) = [linspace(0, 1, floor(kp*5)); linspace(0, 1, floor(kp*5)); linspace(1, 1, floor(kp*5))]';
cmap(floor(kp*5):floor(kp*5)+1,:) = [1 1 1; 1 1 1]; % White color for mid-range
cmap(floor(kp*5)+2:256,:) = [linspace(1, 1, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1)); linspace(1, 0, 256-(floor(kp*5)+1))]';
ax = gca;
ax.Colormap = cmap; % Apply colormap
