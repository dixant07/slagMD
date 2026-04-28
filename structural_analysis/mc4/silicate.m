% =========================================================================
% Silicate Q^n Species Analysis via Adjacency Matrices
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms

% 1. Define Atom Types
type_Si = 2;
type_Al = 4;
type_O  = 5;

% 2. Define Rigorous Cutoff Radii (From your earlier g(r) analysis)
% Do NOT guess these. Use the first minimum of your partial RDFs.
R_cut_SiO = 2.35; 
R_cut_AlO = 2.55; 

% 3. Subsampling Stride
stride = 1000; 

% Trackers for Qn species (Q0, Q1, Q2, Q3, Q4, and 'Other' for anomalies)
Qn_counts = zeros(1, 6); % Indices 1-5 map to Q0-Q4, Index 6 is 'Other'
total_Si_analyzed = 0;
processed_frames = 0;

% =========================================================================
% Main Parsing Loop
% =========================================================================
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file.');
end

while ~feof(fid)
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    
    if contains(tline, 'ITEM: TIMESTEP')
        timestep = str2double(fgetl(fid));
        
        fgetl(fid); 
        N = str2double(fgetl(fid));
        
        % Subsampling check
        if mod(timestep, stride) ~= 0
            for skip = 1:(N+6); fgetl(fid); end
            continue; 
        end
        
        fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); % Skip Box Bounds
        fgetl(fid); % Skip ITEM: ATOMS header
        
        % Read the coordinate block
        data = textscan(fid, '%f %f %f %f %f %f', N);
        types = data{2};
        coords = [data{4}, data{5}, data{6}];
        
        % Isolate Coordinates by Type
        coords_Si = coords(types == type_Si, :);
        coords_Al = coords(types == type_Al, :);
        coords_O  = coords(types == type_O, :);
        
        N_Si = size(coords_Si, 1);
        N_Al = size(coords_Al, 1);
        N_O  = size(coords_O, 1);
        
        % ===============================================================
        % STEP 1: Build Boolean Adjacency Matrices
        % ===============================================================
        % Initialize matrices
        A_SiO = false(N_O, N_Si);
        A_AlO = false(N_O, N_Al);
        
        % Build Si-O Adjacency
        for j = 1:N_Si
            delta = coords_O - coords_Si(j, :);
            delta = delta - L * round(delta / L); % Minimum Image
            r_sq = sum(delta.^2, 2);
            A_SiO(:, j) = r_sq < (R_cut_SiO^2); % Logical array
        end
        
        % Build Al-O Adjacency
        for k = 1:N_Al
            delta = coords_O - coords_Al(k, :);
            delta = delta - L * round(delta / L); % Minimum Image
            r_sq = sum(delta.^2, 2);
            A_AlO(:, k) = r_sq < (R_cut_AlO^2); % Logical array
        end
        
        % ===============================================================
        % STEP 2: Tag Bridging Oxygens (BO)
        % ===============================================================
        % Sum across columns to get total network formers attached to each Oxygen
        Z_O = sum(A_SiO, 2) + sum(A_AlO, 2);
        
        % Create boolean vector: 1 if BO (Z=2), 0 otherwise
        b_BO = (Z_O == 2); 
        
        % ===============================================================
        % STEP 3: The Q^n Dot Product
        % ===============================================================
        % Multiply transposed Si-O adjacency matrix by the BO boolean vector.
        % This instantly yields an N_Si x 1 vector containing the Q^n state (0-4) for every Silicon.
        n_vector = A_SiO' * b_BO; 
        
        % ===============================================================
        % STEP 4: Tally Results for this Frame
        % ===============================================================
        for j = 1:N_Si
            n_val = n_vector(j);
            if n_val >= 0 && n_val <= 4
                Qn_counts(n_val + 1) = Qn_counts(n_val + 1) + 1; % MATLAB is 1-indexed
            else
                Qn_counts(6) = Qn_counts(6) + 1; % 'Other' (Anomalous geometry)
            end
        end
        
        total_Si_analyzed = total_Si_analyzed + N_Si;
        processed_frames = processed_frames + 1;
        
        fprintf('Processed Frame %d (Timestep %d)\n', processed_frames, timestep);
    end
end
fclose(fid);

% =========================================================================
% Final Calculation and Output
% =========================================================================
Qn_percentages = (Qn_counts / total_Si_analyzed) * 100;

fprintf('\n=======================================\n');
fprintf('SILICATE Q^n SPECIES DISTRIBUTION:\n');
fprintf('=======================================\n');
fprintf('Q^0 (Monomer) : %6.2f %%\n', Qn_percentages(1));
fprintf('Q^1 (Dimer)   : %6.2f %%\n', Qn_percentages(2));
fprintf('Q^2 (Chain)   : %6.2f %%\n', Qn_percentages(3));
fprintf('Q^3 (Sheet)   : %6.2f %%\n', Qn_percentages(4));
fprintf('Q^4 (3D Net)  : %6.2f %%\n', Qn_percentages(5));
fprintf('Anomalous (>4): %6.2f %%\n', Qn_percentages(6));
fprintf('=======================================\n');
