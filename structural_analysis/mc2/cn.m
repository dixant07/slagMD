clear all;
clc;
% =========================================================================
% Coordination Number (CN) via Direct Counting from LAMMPS Dump
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms

% 1. Define Target Atom Types (1:Ca, 2:Si, 3:Mg, 4:Al, 5:O)
typeA = 4; % Central atom (e.g., 2 for Silicon, 4 for Aluminum)
typeB = 5; % Neighbor atom (e.g., 5 for Oxygen)

% 2. Define Coordination Cutoff Radius (R_cut)
% YOU MUST GET THIS FROM YOUR RDF PLOT! (The first valley)
% Typical value for Si-O is ~2.35 A, for Al-O is ~2.55 A.
R_cut = 2.55; 

% 3. Subsampling Stride to save execution time
stride = 1000; 

% Trackers for averaging
total_neighbors_found = 0;
total_A_atoms_analyzed = 0;
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
        
        fgetl(fid); % Read number of atoms header
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
        
        % ===============================================================
        % Isolate Coordinates by Type
        % ===============================================================
        coordsA = coords(types == typeA, :);
        coordsB = coords(types == typeB, :);
        
        Na = size(coordsA, 1);
        
        % Tracker for this specific frame
        neighbors_in_this_frame = 0;
        
        % ===============================================================
        % Calculate Distances and Count Neighbors
        % ===============================================================
        for i = 1:Na
            % Distance from atom A(i) to all atoms of type B
            delta = coordsB - coordsA(i, :);
            
            % Minimum Image Convention (Periodic Boundaries)
            delta = delta - L * round(delta / L);
            
            % Scalar distances
            r = sqrt(sum(delta.^2, 2));
            
            % DIRECT COUNTING: How many B atoms are closer than R_cut?
            % (r > 0 ensures we don't count the atom itself if typeA == typeB)
            num_neighbors = sum(r > 0 & r < R_cut);
            
            neighbors_in_this_frame = neighbors_in_this_frame + num_neighbors;
        end
        
        % Accumulate statistics
        total_neighbors_found = total_neighbors_found + neighbors_in_this_frame;
        total_A_atoms_analyzed = total_A_atoms_analyzed + Na;
        processed_frames = processed_frames + 1;
        
        % Optional: Print frame-by-frame CN to monitor convergence
        frame_CN = neighbors_in_this_frame / Na;
        fprintf('Frame %d (Timestep %d): CN = %.4f\n', processed_frames, timestep, frame_CN);
    end
end
fclose(fid);

% =========================================================================
% Final Calculation
% =========================================================================
% Average Coordination Number across all analyzed atoms and frames
final_CN = total_neighbors_found / total_A_atoms_analyzed;

fprintf('\n=======================================\n');
fprintf('FINAL COORDINATION NUMBER RESULTS:\n');
fprintf('Central Atom Type : %d\n', typeA);
fprintf('Neighbor Atom Type: %d\n', typeB);
fprintf('Cutoff Radius     : %.2f A\n', R_cut);
fprintf('Average CN        : %.4f\n', final_CN);
fprintf('=======================================\n');
