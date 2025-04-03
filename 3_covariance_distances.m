% Generate TCL file extractor file

% Set variables
% Define indices for RBD and NB regions
rbd_ind = 334:528;
nb_ind = 1:128;

% Define interaction interface indices for RBD and NB
rbd_int = 445:506;
nb_int = [26:32 52:56 99:116];

% Calculate the size of RBD and NB regions
rbd_size = size(rbd_ind,2);
nb_size = size(nb_ind,2);

% Define residue types
hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "CYS", "PRO"];
acidic = ["ASP", "GLU"];
basic = ["LYS", "ARG"];
polar = ["SER", "THR", "ASN", "GLN",  "TYR"];

% Load PSF files for different variants
wt_psf = read_psf('Data/WT/wt_aligned.psf');
alpha_psf = read_psf('Data/Alpha/alpha_aligned.psf');
beta_psf = read_psf('Data/Beta/beta_aligned.psf');
omicron_psf = read_psf('Data/Omicron/omicron_aligned.psf');

% Find attractive pairs for WT variant
[i1 n1] = find(cov_wt_d > 0);

% Write TCL file for WT attractive pairs
fid = fopen('wt_attractive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

% Loop through attractive pairs and write commands to TCL file
for k = 1:length(i1)
    res1 = wt_psf.residue_name{i1(k)};
    res2 = wt_psf.residue_name{n1(k) + rbd_size};
    
    % Check residue types and write appropriate commands
    if ismember(res1, hydrophobic) && ismember(res2, hydrophobic)
        fprintf(fid, 'set s1 "sidechain and noh and segname %s and resid %d"\n', wt_psf.segment_name{i1(k)}, wt_psf.residue_id(i1(k)));
        fprintf(fid, 'set s2 "sidechain and noh and segname %s and resid %d"\n', wt_psf.segment_name{n1(k) + rbd_size}, wt_psf.residue_id(n1(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, acidic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', wt_psf.segment_name{i1(k)}, wt_psf.residue_id(i1(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', wt_psf.segment_name{n1(k) + rbd_size}, wt_psf.residue_id(n1(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', wt_psf.segment_name{i1(k)}, wt_psf.residue_id(i1(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', wt_psf.segment_name{n1(k) + rbd_size}, wt_psf.residue_id(n1(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
end
fclose(fid);

% Repeat the same process for other variants (Alpha, Beta, Omicron)
% Find attractive pairs for Alpha variant
[i3 n3] = find(cov_alpha_d > 0);

% Write TCL file for Alpha attractive pairs
fid = fopen('alpha_attractive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i3)
    res1 = alpha_psf.residue_name{i3(k)};
    res2 = alpha_psf.residue_name{n3(k) + rbd_size};
    
    if ismember(res1, hydrophobic) && ismember(res2, hydrophobic)
        fprintf(fid, 'set s1 "sidechain and noh and segname %s and resid %d"\n', alpha_psf.segment_name{i3(k)}, alpha_psf.residue_id(i3(k)));
        fprintf(fid, 'set s2 "sidechain and noh and segname %s and resid %d"\n', alpha_psf.segment_name{n3(k) + rbd_size}, alpha_psf.residue_id(n3(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, acidic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', alpha_psf.segment_name{i3(k)}, alpha_psf.residue_id(i3(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', alpha_psf.segment_name{n3(k) + rbd_size}, alpha_psf.residue_id(n3(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', alpha_psf.segment_name{i3(k)}, alpha_psf.residue_id(i3(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', alpha_psf.segment_name{n3(k) + rbd_size}, alpha_psf.residue_id(n3(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
end
fclose(fid);

% Repeat similar logic for Beta and Omicron variants for both attractive and repulsive pairs
% The rest of the code follows the same structure for other variants and pair types
