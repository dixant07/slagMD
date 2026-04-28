clear all
clc

% =========================================================================
% Partial Radial Distribution Function g_AB(r) with Subsampling
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms
rmax = 22.0;                  % Maximum radius (L/2)
dr = 0.1;                     % Bin size in Angstroms

% 1. Target Atom Types (Based on: 1:Ca, 2:Si, 3:Mg, 4:Al, 5:O)
typeA = 2; % Central atom (e.g., Silicon)
typeB = 5; % Target atom (e.g., Oxygen)

% 2. Statistical Subsampling (The Speed Fix)
% For 300,000 frames, a stride of 1000 gives 300 independent snapshots.
stride = 1000; 

% Setup Histogram
nbins = floor(rmax / dr);
g_hist = zeros(nbins, 1);
r_axis = (0.5:1:nbins)' * dr;
processed_frames = 0;

% Track the number of atoms of type A and B for normalization
avg_Na = 0;
avg_Nb = 0;

% 3. Open File and Parse
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file.');
end

while ~feof(fid)
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    
    if contains(tline, 'ITEM: TIMESTEP')
        timestep = str2double(fgetl(fid));
        
        % Read number of atoms
        fgetl(fid); 
        N = str2double(fgetl(fid));
        
        % Check if we should process or skip this frame
        if mod(timestep, stride) ~= 0
            % Skip the remaining (N + 6) lines of this frame to save massive time
            for skip = 1:(N+6)
                fgetl(fid);
            end
            continue; 
        end
        
        % If we reach here, we are processing the frame
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
        Nb = size(coordsB, 1);
        
        avg_Na = avg_Na + Na;
        avg_Nb = avg_Nb + Nb;
        
        % ===============================================================
        % Calculate Distances (Cross-Pair Logic)
        % ===============================================================
        for i = 1:Na
            if typeA == typeB
                % If computing A-A (e.g., O-O), avoid double counting
                delta = coordsB((i+1):end, :) - coordsA(i, :);
                add_val = 2; % Add 2 because i->j and j->i are found simultaneously
            else
                % If computing A-B (e.g., Si-O), calculate all A to all B
                delta = coordsB - coordsA(i, :);
                add_val = 1; % Add 1 because we are strictly looking from A to B
            end
            
            % Minimum Image Convention
            delta = delta - L * round(delta / L);
            
            % Scalar distances
            r = sqrt(sum(delta.^2, 2));
            
            % Bin valid distances
            r_valid = r(r < rmax);
            bin_indices = floor(r_valid / dr) + 1;
            
            for b = 1:length(bin_indices)
                idx = bin_indices(b);
                if idx <= nbins
                    g_hist(idx) = g_hist(idx) + add_val;
                end
            end
        end
        
        processed_frames = processed_frames + 1;
        fprintf('Processed sampled frame: %d (Timestep %d)\n', processed_frames, timestep);
    end
end
fclose(fid);