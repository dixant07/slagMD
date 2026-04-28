% =========================================================================
% Radial Distribution Function g(r) from LAMMPS Dump
% =========================================================================

% 1. System Parameters
filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms
rmax = 22.0;                  % Maximum radius (L/2)
dr = 0.1;                     % Bin size in Angstroms

% 2. Setup Histogram
nbins = floor(rmax / dr);
g_hist = zeros(nbins, 1);
r_axis = (0.5:1:nbins)' * dr; % Center of each bin
total_frames = 0;

% 3. Open File and Parse
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file.');
end

while ~feof(fid)
    % Read header lines
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    
    if contains(tline, 'ITEM: TIMESTEP')
        timestep = str2double(fgetl(fid));
        
        % Read number of atoms
        fgetl(fid); % ITEM: NUMBER OF ATOMS
        N = str2double(fgetl(fid));
        
        % Skip Box Bounds (We already defined L = 44)
        fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
        
        % Read ITEM: ATOMS id type q x y z
        fgetl(fid); 
        
        % Read the coordinate block using textscan for speed
        % Format: %f %f %f %f %f %f (id, type, q, x, y, z)
        data = textscan(fid, '%f %f %f %f %f %f', N);
        
        % Extract X, Y, Z coordinates into a single Nx3 matrix
        coords = [data{4}, data{5}, data{6}];
        
        % ===============================================================
        % Calculate Distances (Vectorized for MATLAB Speed)
        % ===============================================================
        for i = 1:(N-1)
            % Distance from atom i to ALL subsequent atoms j
            % This array subtraction replaces the inner loop
            delta = coords((i+1):end, :) - coords(i, :);
            
            % Apply Minimum Image Convention
            delta = delta - L * round(delta / L);
            
            % Calculate scalar distances (r = sqrt(dx^2 + dy^2 + dz^2))
            r = sqrt(sum(delta.^2, 2));
            
            % Keep only pairs within our cutoff
            r_valid = r(r < rmax);
            
            % Bin the valid distances
            % floor(r / dr) + 1 gives the correct bin index in MATLAB (1-based)
            bin_indices = floor(r_valid / dr) + 1;
            
            % Populate the histogram (adding 2 because pair counts for both i and j)
            for b = 1:length(bin_indices)
                idx = bin_indices(b);
                if idx <= nbins
                    g_hist(idx) = g_hist(idx) + 2;
                end
            end
        end
        
        total_frames = total_frames + 1;
        fprintf('Processed frame: %d\n', total_frames);
    end
end
fclose(fid);

% =========================================================================
% 4. Normalization to ideal gas density
% =========================================================================
rho = N / (L^3); % Global number density

g_r = zeros(nbins, 1);
for k = 1:nbins
    r_inner = (k - 1) * dr;
    r_outer = k * dr;
    
    % Exact volume of the spherical shell
    V_shell = (4/3) * pi * (r_outer^3 - r_inner^3);
    
    % Ideal number of atoms in this shell
    n_ideal = rho * V_shell;
    
    % Normalize histogram
    % g(r) = count / (frames * N * n_ideal)
    g_r(k) = g_hist(k) / (total_frames * N * n_ideal);
end

% =========================================================================
% 5. Plot the result
% =========================================================================
figure;
plot(r_axis, g_r, 'b-', 'LineWidth', 1.5);
xlabel('Distance r (Angstroms)', 'FontSize', 12);
ylabel('g(r)', 'FontSize', 12);
title('Total Radial Distribution Function', 'FontSize', 14);
grid on;
xlim([0 rmax]);
