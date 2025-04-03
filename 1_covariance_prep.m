% Title: DCD Data Reader and Covariance Matrix Constructor

% Authors: Mert Golcuk, PhD (Postdoctoral Researcher), Mert Gur, PhD (PI)

% Lab: Gur Lab, University of Pittsburgh School of Medicine
% Date: 2025

% Description:
% This script reads molecular dynamics trajectory data from DCD files and constructs a covariance matrix
% based on atomic positional fluctuations. The resulting covariance matrix is typically used for
% principal component analysis (PCA), structural dynamics analysis, and related computational biophysics studies.

% Requirements:
% - matdcd package

% Usage:
% - Specify the input DCD file.
% - Run the script to compute and output the covariance matrix.

% License:
% Available for academic and research purposes under MIT License.


clear all; close all; clc

% Add the matdcd package to the MATLAB path
addpath("matdcd-1.0/")

% Set cutoff values for attraction and repulsion
coff_att = 11; % Attraction cutoff
coff_rep = 14; % Repulsion cutoff

% Define residue indices for RBD and N-terminal domain
rbd_ind = 334:528; % RBD indices
nb_ind = 1:128; % N-terminal domain indices

% Define interaction regions for RBD and N-terminal domain
rbd_int = 445:506; % RBD interaction region
nb_int = [26:32 52:56 99:116]; % N-terminal domain interaction region

% Load DCD files for different variants
wt = readdcd("Data/wt_aligned.dcd",1:323); % Wild-type
alpha = readdcd("Data/n501y_aligned.dcd",1:323); % Alpha variant
beta = readdcd("Data/triple_aligned.dcd",1:323); % Beta variant
omicron = readdcd("Data/omicron_aligned.dcd",1:323); % Omicron variant

% Compute covariance matrices for each variant
cov_wt = compute_covariance(wt);
cov_alpha = compute_covariance(alpha);
cov_beta = compute_covariance(beta);
cov_omicron = compute_covariance(omicron);

% Extract cross-covariance matrices between RBD and N-terminal domain
cov_wt_a = cov_wt(1:size(rbd_ind,2),size(rbd_ind,2)+1:end);
cov_alpha_a = cov_alpha(1:size(rbd_ind,2),size(rbd_ind,2)+1:end);
cov_beta_a = cov_beta(1:size(rbd_ind,2),size(rbd_ind,2)+1:end);
cov_omicron_a = cov_omicron(1:size(rbd_ind,2),size(rbd_ind,2)+1:end);

% Initialize matrices for further processing
cov_wt_b = cov_wt_a;
cov_alpha_b = cov_alpha_a;
cov_beta_b = cov_beta_a;
cov_omicron_b = cov_omicron_a;

cov_wt_c = cov_wt_a;
cov_alpha_c = cov_alpha_a;
cov_beta_c = cov_beta_a;
cov_omicron_c = cov_omicron_a;

% Initialize matrices to store final results
cov_wt_d = zeros(size(cov_wt_a,1), size(cov_wt_a,2));
cov_alpha_d = zeros(size(cov_alpha_a,1), size(cov_alpha_a,2));
cov_beta_d = zeros(size(cov_beta_a,1), size(cov_beta_a,2));
cov_omicron_d = zeros(size(cov_omicron_a,1), size(cov_omicron_a,2));

% Compute distance matrices for attraction cutoff
wt_dist = dist_mat_xyz(mean(wt),coff_att);
alpha_dist = dist_mat_xyz(mean(alpha),coff_att);
beta_dist = dist_mat_xyz(mean(beta),coff_att);
omicron_dist = dist_mat_xyz(mean(omicron),coff_att);

% Compute distance matrices for repulsion cutoff
wt_dist2 = dist_mat_xyz(mean(wt),coff_rep);
alpha_dist2 = dist_mat_xyz(mean(alpha),coff_rep);
beta_dist2 = dist_mat_xyz(mean(beta),coff_rep);
omicron_dist2 = dist_mat_xyz(mean(omicron),coff_rep);

% Set covariance values to zero where distances exceed attraction cutoff
cov_wt_b(wt_dist(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_alpha_b(alpha_dist(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_beta_b(beta_dist(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_omicron_b(omicron_dist(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;

% Set covariance values to zero where distances exceed repulsion cutoff
cov_wt_c(wt_dist2(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_alpha_c(alpha_dist2(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_beta_c(beta_dist2(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;
cov_omicron_c(omicron_dist2(1:size(rbd_ind,2),size(rbd_ind,2)+1:end) == 0) = 0;

% Identify indices where covariance values meet specific conditions
[i1 n1] = find(cov_wt_c < 0); % Negative values in cov_wt_c
[i2 n2] = find(cov_wt_b > 0); % Positive values in cov_wt_b
[i3 n3] = find(cov_alpha_c < 0); % Negative values in cov_alpha_c
[i4 n4] = find(cov_alpha_b > 0); % Positive values in cov_alpha_b
[i5 n5] = find(cov_beta_c < 0); % Negative values in cov_beta_c
[i6 n6] = find(cov_beta_b > 0); % Positive values in cov_beta_b
[i7 n7] = find(cov_omicron_c < 0); % Negative values in cov_omicron_c
[i8 n8] = find(cov_omicron_b > 0); % Positive values in cov_omicron_b

% Update final covariance matrices based on identified indices
cov_wt_d(i1,n1) = cov_wt_c(i1,n1);
cov_wt_d(i2,n2) = cov_wt_b(i2,n2);
cov_alpha_d(i3,n3) = cov_alpha_c(i3,n3);
cov_alpha_d(i4,n4) = cov_alpha_b(i4,n4);
cov_beta_d(i5,n5) = cov_beta_c(i5,n5);
cov_beta_d(i6,n6) = cov_beta_b(i6,n6);
cov_omicron_d(i7,n7) = cov_omicron_c(i7,n7);
cov_omicron_d(i8,n8) = cov_omicron_b(i8,n8);
