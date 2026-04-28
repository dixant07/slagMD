% =========================================================================
% Oxygen Network Topology Analysis (BO, NBO, FO)
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms

% 1. Define Atom Types
type_Si = 2;
type_Al = 4;
type_O  = 5;

% 2. Define Cutoff Radii (From your RDF first-valley minimums)
% Adjust these slightly based on your exact g(r) plots if needed.
R_cut_SiO = 2.35; 
R_cut_AlO = 2.55; 

% 3. Subsampling Stride
stride = 1000; 

% Trackers for averaging percentages across frames
sum_pct_FO  = 0;
sum_pct_NBO = 0;
sum_pct_BO  = 0;
sum_pct_Tri = 0;
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
        
        num_O = size(coords_O, 1);
        
        % Frame-specific counters
        count_FO  = 0;
        count_NBO = 0;
        count_BO  = 0;
        count_Tri = 0;
        
        % ===============================================================
        % Calculate Bonding for Every Oxygen Atom
        % ===============================================================
        for i = 1:num_O
            
            % 1. Count bonds to Silicon
            delta_Si = coords_Si - coords_O(i, :);
            delta_Si = delta_Si - L * round(delta_Si / L); % Minimum Image
            r_Si = sqrt(sum(delta_Si.^2, 2));
            bonds_Si = sum(r_Si < R_cut_SiO);
            
            % 2. Count bonds to Aluminum
            delta_Al = coords_Al - coords_O(i, :);
            delta_Al = delta_Al - L * round(delta_Al / L); % Minimum Image
            r_Al = sqrt(sum(delta_Al.^2, 2));
            bonds_Al = sum(r_Al < R_cut_AlO);
            
            % 3. Total Network Former Coordination (Z)
            Z = bonds_Si + bonds_Al;
            
            % 4. Classify the Oxygen
            if Z == 0
                count_FO = count_FO + 1;
            elseif Z == 1
                count_NBO = count_NBO + 1;
            elseif Z == 2
                count_BO = count_BO + 1;
            else
                count_Tri = count_Tri + 1; % Z >= 3 (Triclusters)
            end
        end
        
        % Calculate percentages for this frame
        pct_FO  = (count_FO / num_O) * 100;
        pct_NBO = (count_NBO / num_O) * 100;
        pct_BO  = (count_BO / num_O) * 100;
        pct_Tri = (count_Tri / num_O) * 100;
        
        % Accumulate for final averaging
        sum_pct_FO  = sum_pct_FO  + pct_FO;
        sum_pct_NBO = sum_pct_NBO + pct_NBO;
        sum_pct_BO  = sum_pct_BO  + pct_BO;
        sum_pct_Tri = sum_pct_Tri + pct_Tri;
        
        processed_frames = processed_frames + 1;
        
        fprintf('Frame %d | BO: %.2f%%, NBO: %.2f%%, FO: %.2f%%\n', ...
                processed_frames, pct_BO, pct_NBO, pct_FO);
    end
end
fclose(fid);

% =========================================================================
% Final Calculation
% =========================================================================
final_FO  = sum_pct_FO / processed_frames;
final_NBO = sum_pct_NBO / processed_frames;
final_BO  = sum_pct_BO / processed_frames;
final_Tri = sum_pct_Tri / processed_frames;

fprintf('\n=======================================\n');
fprintf('OXYGEN NETWORK TOPOLOGY (Averaged):\n');
fprintf('=======================================\n');
fprintf('Free Oxygen (FO)         [Z=0]: %6.2f %%\n', final_FO);
fprintf('Non-Bridge Oxygen (NBO)  [Z=1]: %6.2f %%\n', final_NBO);
fprintf('Bridge Oxygen (BO)       [Z=2]: %6.2f %%\n', final_BO);
fprintf('Oxygen Triclusters       [Z>2]: %6.2f %%\n', final_Tri);
fprintf('=======================================\n');
