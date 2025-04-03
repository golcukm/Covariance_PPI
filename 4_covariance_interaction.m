%close all

% load WT interaction data

dataFiles = dir('Data/WT/distance*.dat');
wt_filenames = cell(length(dataFiles), 1);
numFiles = length(dataFiles);

% Load the first file to determine the number of rows
firstFileName = fullfile(dataFiles(1).folder, dataFiles(1).name);
firstData = load(firstFileName);
numRows = size(firstData, 1);

% Preallocate the matrix for second columns
wt_int = zeros(numRows, numFiles);

for k = 1:numFiles
    fileName = fullfile(dataFiles(k).folder, dataFiles(k).name);
    data = load(fileName);
    wt_filenames{k} = dataFiles(k).name;
    wt_int(:, k) = data(:, 2);
end

% Now fileNames contains the names of the files
% and secondColumns contains the second columns of the data in a matrix

% load Alpha interaction data

dataFiles = dir('Data/N501Y/distance*.dat');
alpha_filenames = cell(length(dataFiles), 1);
numFiles = length(dataFiles);

% Load the first file to determine the number of rows
firstFileName = fullfile(dataFiles(1).folder, dataFiles(1).name);
firstData = load(firstFileName);
numRows = size(firstData, 1);

% Preallocate the matrix for second columns
alpha_int = zeros(numRows, numFiles);

for k = 1:numFiles
    fileName = fullfile(dataFiles(k).folder, dataFiles(k).name);
    data = load(fileName);
    alpha_filenames{k} = dataFiles(k).name;
    alpha_int(:, k) = data(:, 2);
end

% Now fileNames contains the names of the files
% and secondColumns contains the second columns of the data in a matrix

% load Beta interaction data

dataFiles = dir('Data/Triple/distance*.dat');
beta_filenames = cell(length(dataFiles), 1);
numFiles = length(dataFiles);

% Load the first file to determine the number of rows
firstFileName = fullfile(dataFiles(1).folder, dataFiles(1).name);
firstData = load(firstFileName);
numRows = size(firstData, 1);

% Preallocate the matrix for second columns
beta_int = zeros(numRows, numFiles);

for k = 1:numFiles
    fileName = fullfile(dataFiles(k).folder, dataFiles(k).name);
    data = load(fileName);
    beta_filenames{k} = dataFiles(k).name;
    beta_int(:, k) = data(:, 2);
end

% Now fileNames contains the names of the files
% and secondColumns contains the second columns of the data in a matrix

% load Omicron interaction data

dataFiles = dir('Data/Omicron/distance*.dat');
omicron_filenames = cell(length(dataFiles), 1);
numFiles = length(dataFiles);

% Load the first file to determine the number of rows
firstFileName = fullfile(dataFiles(1).folder, dataFiles(1).name);
firstData = load(firstFileName);
numRows = size(firstData, 1);

% Preallocate the matrix for second columns
omicron_int = zeros(numRows, numFiles);

for k = 1:numFiles
    fileName = fullfile(dataFiles(k).folder, dataFiles(k).name);
    data = load(fileName);
    omicron_filenames{k} = dataFiles(k).name;
    omicron_int(:, k) = data(:, 2);
end

% Now fileNames contains the names of the files
% and secondColumns contains the second columns of the data in a matrix


wt_hyp = [1:4 6:size(wt_int, 2)];
alpha_hyp = [1:3 5:size(alpha_int,2)];
beta_hyp = [1:size(beta_int,2)];
omicron_hyp = [1:size(omicron_int,2)];

wt_sb = 5;
alphe_wt = 4;

wt_int_2 = wt_int;
wt_int_2(:, wt_hyp) = wt_int(:, wt_hyp)<=8;
wt_int_2(:, wt_sb) = wt_int(:, wt_sb)<=6;

wt_perc = 100*sum(wt_int_2, 1)/size(wt_int_2, 1);

alpha_int_2 = alpha_int;
alpha_int_2(:, alpha_hyp) = alpha_int(:, alpha_hyp)<=8;
alpha_int_2(:, alphe_wt) = alpha_int(:, alphe_wt)<=6;

alpha_perc = 100*sum(alpha_int_2, 1)/size(alpha_int_2, 1);

beta_int_2 = beta_int;
beta_int_2(:, beta_hyp) = beta_int(:, beta_hyp)<=8;

beta_perc = 100*sum(beta_int_2, 1)/size(beta_int_2, 1);

omicron_int_2 = omicron_int;
omicron_int_2(:, omicron_hyp) = omicron_int(:, omicron_hyp)<=8;

omicron_perc = 100*sum(omicron_int_2, 1)/size(omicron_int_2, 1);


% Extract residue number pairs from wt_filenames
wt_residue_pairs = zeros(length(wt_filenames), 2);

for i = 1:length(wt_filenames)
    % Extract the numbers from the filename
    tokens = regexp(wt_filenames{i}, '_(\d+)_C_.*_(\d+)\.dat', 'tokens');
    if ~isempty(tokens)
        wt_residue_pairs(i, :) = str2double(tokens{1});
    end
end

% Extract residue number pairs from alpha_filenames
alpha_residue_pairs = zeros(length(alpha_filenames), 2);

for i = 1:length(alpha_filenames)
    % Extract the numbers from the filename
    tokens = regexp(alpha_filenames{i}, '_(\d+)_C_.*_(\d+)\.dat', 'tokens');
    if ~isempty(tokens)
        alpha_residue_pairs(i, :) = str2double(tokens{1});
    end
end

% Extract residue number pairs from beta_filenames
beta_residue_pairs = zeros(length(beta_filenames), 2);

for i = 1:length(beta_filenames)
    % Extract the numbers from the filename
    tokens = regexp(beta_filenames{i}, '_(\d+)_C_.*_(\d+)\.dat', 'tokens');
    if ~isempty(tokens)
        beta_residue_pairs(i, :) = str2double(tokens{1});
    end
end

% Extract residue number pairs from omicron_filenames
omicron_residue_pairs = zeros(length(omicron_filenames), 2);

for i = 1:length(omicron_filenames)
    % Extract the numbers from the filename
    tokens = regexp(omicron_filenames{i}, '_(\d+)_C_.*_(\d+)\.dat', 'tokens');
    if ~isempty(tokens)
        omicron_residue_pairs(i, :) = str2double(tokens{1});
    end
end


wt_dist_cov = zeros(size(rbd_ind,2),(size(nb_ind,2)));

for k = 1:size(wt_perc,2)
    wt_dist_cov(wt_residue_pairs(k,1)-333,wt_residue_pairs(k,2)) = wt_perc(k);
end

alpha_dist_cov = zeros(size(rbd_ind,2),(size(nb_ind,2)));  

for k = 1:size(alpha_perc,2)
    alpha_dist_cov(alpha_residue_pairs(k,1)-333,alpha_residue_pairs(k,2)) = alpha_perc(k);
end

beta_dist_cov = zeros(size(rbd_ind,2),(size(nb_ind,2)));
beta_perc(15)=-beta_perc(15);
for k = 1:size(beta_perc,2)
    beta_dist_cov(beta_residue_pairs(k,1)-333,beta_residue_pairs(k,2)) = beta_perc(k);
end

omicron_dist_cov = zeros(size(rbd_ind,2),(size(nb_ind,2)));

for k = 1:size(omicron_perc,2)
    omicron_dist_cov(omicron_residue_pairs(k,1)-333,omicron_residue_pairs(k,2)) = omicron_perc(k);
end
