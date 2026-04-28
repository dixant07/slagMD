% =========================================================================
% Aluminum Coordination Analysis ([AlO4], [AlO5])
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Your specific trajectory file
L = 44.0;                         % Box size in Angstroms

% 1. Define Atom Types (1:Ca, 2:Si, 3:Mg, 4:Al, 5:O)
type_Al = 4;
type_O  = 5;

% 2. Rigorous Cutoff Radii
R_cut_SiO = 2.35; % Maintained for your records/consistency
R_cut_AlO = 2.55; % Used for this calculation

% 3. Subsampling Stride (Adjust if you want higher/lower temporal resolution)
stride = 1000; 

% Trackers for Aluminum Coordination
count_AlO4 = 0;
count_AlO5 = 0;
count_AlO_Other = 0; % Captures 3-coordinated, 6-coordinated, etc.
total_Al_analyzed = 0;
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
        
        % Read the coordinate block using fast textscan
        data = textscan(fid, '%f %f %f %f %f %f', N);
        types = data{2};
        coords = [data{4}, data{5}, data{6}];
        
        % Isolate Coordinates by Type
        coords_Al = coords(types == type_Al, :);
        coords_O  = coords(types == type_O, :);
        
        N_Al = size(coords_Al, 1);
        
        % Frame-specific counters
        frame_AlO4 = 0;
        frame_AlO5 = 0;
        frame_Other = 0;
        
        % ===============================================================
        % Calculate Coordination for Every Aluminum Atom
        % ===============================================================
        % We use the squared cutoff to avoid computing expensive square roots
        R_cut_AlO_sq = R_cut_AlO^2; 
        
        for k = 1:N_Al
            % Vectorized distance from Aluminum k to ALL Oxygen atoms
            delta = coords_O - coords_Al(k, :);
            
            % Minimum Image Convention (Periodic Boundaries)
            delta = delta - L * round(delta / L);
            
            % Squared scalar distance
            r_sq = sum(delta.^2, 2);
            
            % Count how many Oxygens are within the cutoff
            Z_Al = sum(r_sq < R_cut_AlO_sq);
            
            % Tally the coordination state
            if Z_Al == 4
                frame_AlO4 = frame_AlO4 + 1;
            elseif Z_Al == 5
                frame_AlO5 = frame_AlO5 + 1;
            else
                frame_Other = frame_Other + 1;
            end
        end
        
        % Accumulate statistics across frames
        count_AlO4 = count_AlO4 + frame_AlO4;
        count_AlO5 = count_AlO5 + frame_AlO5;
        count_AlO_Other = count_AlO_Other + frame_Other;
        total_Al_analyzed = total_Al_analyzed + N_Al;
        
        processed_frames = processed_frames + 1;
        
        fprintf('Processed Frame %d (Timestep %d) | [AlO4]: %d, [AlO5]: %d\n', ...
                processed_frames, timestep, frame_AlO4, frame_AlO5);
    end
end
fclose(fid);

% =========================================================================
% Final Calculation and Output
% =========================================================================
pct_AlO4  = (count_AlO4 / total_Al_analyzed) * 100;
pct_AlO5  = (count_AlO5 / total_Al_analyzed) * 100;
pct_Other = (count_AlO_Other / total_Al_analyzed) * 100;

fprintf('\n=======================================\n');
fprintf('ALUMINUM COORDINATION DISTRIBUTION:\n');
fprintf('=======================================\n');
fprintf('[AlO4] (Tetrahedral)  : %6.2f %%\n', pct_AlO4);
fprintf('[AlO5] (Pentahedral)  : %6.2f %%\n', pct_AlO5);
fprintf('Other  (3- or 6-coord): %6.2f %%\n', pct_Other);
fprintf('=======================================\n');
