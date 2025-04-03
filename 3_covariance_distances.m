% generate tcl file extractor file

% Set variables
% RBD and NB indices
rbd_ind = 334:528;
nb_ind = 1:128;

% RBD and NB interaction interface indices
rbd_int = 445:506;
nb_int = [26:32 52:56 99:116];

% Size of RBD and NB
rbd_size = size(rbd_ind,2);
nb_size = size(nb_ind,2);

% Set residue types
hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "CYS", "PRO"];
acidic = ["ASP", "GLU"];
basic = ["LYS", "ARG"];
polar = ["SER", "THR", "ASN", "GLN",  "TYR"];

% Laod psf files
wt_psf = read_psf('Data/wt_aligned.psf');
alpha_psf = read_psf('Data/n501y_aligned.psf');
beta_psf = read_psf('Data/triple_aligned.psf');
omicron_psf = read_psf('Data/omicron_aligned.psf');

% Find attractive pairs
[i1 n1] = find(cov_wt_d > 0);

% Write tcl file
fid = fopen('wt_attractive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i1)
    res1 = wt_psf.residue_name{i1(k)};
    res2 = wt_psf.residue_name{n1(k) + rbd_size};
    
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

[i3 n3] = find(cov_alpha_d > 0);

% Write tcl file
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

[i5 n5] = find(cov_beta_d > 0);

% Write tcl file
fid = fopen('beta_attractive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i5)
    res1 = beta_psf.residue_name{i5(k)};
    res2 = beta_psf.residue_name{n5(k) + rbd_size};
    
    if ismember(res1, hydrophobic) && ismember(res2, hydrophobic)
        fprintf(fid, 'set s1 "sidechain and noh and segname %s and resid %d"\n', beta_psf.segment_name{i5(k)}, beta_psf.residue_id(i5(k)));
        fprintf(fid, 'set s2 "sidechain and noh and segname %s and resid %d"\n', beta_psf.segment_name{n5(k) + rbd_size}, beta_psf.residue_id(n5(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, acidic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', beta_psf.segment_name{i5(k)}, beta_psf.residue_id(i5(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', beta_psf.segment_name{n5(k) + rbd_size}, beta_psf.residue_id(n5(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', beta_psf.segment_name{i5(k)}, beta_psf.residue_id(i5(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', beta_psf.segment_name{n5(k) + rbd_size}, beta_psf.residue_id(n5(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);

[i7 n7] = find(cov_omicron_d > 0);

% Write tcl file
fid = fopen('omicron_attractive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i7)
    res1 = omicron_psf.residue_name{i7(k)};
    res2 = omicron_psf.residue_name{n7(k) + rbd_size};
    
    if ismember(res1, hydrophobic) && ismember(res2, hydrophobic)
        fprintf(fid, 'set s1 "sidechain and noh and segname %s and resid %d"\n', omicron_psf.segment_name{i7(k)}, omicron_psf.residue_id(i7(k)));
        fprintf(fid, 'set s2 "sidechain and noh and segname %s and resid %d"\n', omicron_psf.segment_name{n7(k) + rbd_size}, omicron_psf.residue_id(n7(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, acidic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', omicron_psf.segment_name{i7(k)}, omicron_psf.residue_id(i7(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', omicron_psf.segment_name{n7(k) + rbd_size}, omicron_psf.residue_id(n7(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', omicron_psf.segment_name{i7(k)}, omicron_psf.residue_id(i7(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', omicron_psf.segment_name{n7(k) + rbd_size}, omicron_psf.residue_id(n7(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);


% Find Repulsive pairs
[i2 n2] = find(cov_wt_d < 0);

% Write tcl file
fid = fopen('wt_repulsive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i2)
    res1 = wt_psf.residue_name{i2(k)};
    res2 = wt_psf.residue_name{n2(k) + rbd_size};
    
    if ismember(res1, acidic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', wt_psf.segment_name{i2(k)}, wt_psf.residue_id(i1(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', wt_psf.segment_name{n2(k) + rbd_size}, wt_psf.residue_id(n2(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', wt_psf.segment_name{i2(k)}, wt_psf.residue_id(i1(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', wt_psf.segment_name{n2(k) + rbd_size}, wt_psf.residue_id(n2(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);

[i4 n4] = find(cov_alpha_d < 0);


% Write tcl file
fid = fopen('alpha_repulsive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i4)
    res1 = alpha_psf.residue_name{i4(k)};
    res2 = alpha_psf.residue_name{n4(k) + rbd_size};
    
    if ismember(res1, acidic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', alpha_psf.segment_name{i4(k)}, alpha_psf.residue_id(i4(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', alpha_psf.segment_name{n4(k) + rbd_size}, alpha_psf.residue_id(n4(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', alpha_psf.segment_name{i4(k)}, alpha_psf.residue_id(i4(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', alpha_psf.segment_name{n4(k) + rbd_size}, alpha_psf.residue_id(n4(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);

[i6 n6] = find(cov_beta_d < 0);

% Write tcl file
fid = fopen('beta_repulsive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i6)
    res1 = beta_psf.residue_name{i6(k)};
    res2 = beta_psf.residue_name{n6(k) + rbd_size};
    
    if ismember(res1, acidic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', beta_psf.segment_name{i6(k)}, beta_psf.residue_id(i6(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', beta_psf.segment_name{n6(k) + rbd_size}, beta_psf.residue_id(n6(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', beta_psf.segment_name{i6(k)}, beta_psf.residue_id(i6(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', beta_psf.segment_name{n6(k) + rbd_size}, beta_psf.residue_id(n6(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);

[i8 n8] = find(cov_omicron_d < 0);

% Write tcl file
fid = fopen('omicron_repulsive.tcl', 'w');
fprintf(fid, 'source distance_measure.tcl\n');

for k = 1:length(i8)
    res1 = omicron_psf.residue_name{i8(k)};
    res2 = omicron_psf.residue_name{n8(k) + rbd_size};
    
    if ismember(res1, acidic) && ismember(res2, acidic)
        fprintf(fid, 'set s1 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', omicron_psf.segment_name{i8(k)}, omicron_psf.residue_id(i8(k)));
        fprintf(fid, 'set s2 "name OE1 OE2 OD1 OD2 and noh and segname %s and resid %d"\n', omicron_psf.segment_name{n8(k) + rbd_size}, omicron_psf.residue_id(n8(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    elseif ismember(res1, basic) && ismember(res2, basic)
        fprintf(fid, 'set s1 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', omicron_psf.segment_name{i8(k)}, omicron_psf.residue_id(i8(k)));
        fprintf(fid, 'set s2 "name NH1 NH2 NZ and noh and segname %s and resid %d"\n', omicron_psf.segment_name{n8(k) + rbd_size}, omicron_psf.residue_id(n8(k) + rbd_size));
        fprintf(fid, 'calculate_distance $s1 $s2\n');
    end
    
end
fclose(fid);
