% =========================================================================
% Total System Mean Square Displacement (MSD)
% =========================================================================

filename = 'dump.slag.lammpstrj'; % Replace with your actual file name
L = 44.0;                     % Box size in Angstroms
dt_fs = 1.0;                  % MD integration timestep (e.g., 1 fs)

% Subsampling Stride
stride = 1000; 

% Initialization flags
is_first_frame = true;
time_array = [];
msd_array = [];

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
        
        % ---------------------------------------------------------------
        % CRITICAL STEP: Sort ALL atoms by ID 
        % ---------------------------------------------------------------
        [~, sort_idx] = sort(data{1}); 
        % We don't need types anymore for total MSD
        coords_all = [data{4}(sort_idx), data{5}(sort_idx), data{6}(sort_idx)];
        
        % ---------------------------------------------------------------
        % Initialize or Unwrap
        % ---------------------------------------------------------------
        if is_first_frame
            r0 = coords_all;
            r_prev = coords_all;
            r_unwrap = coords_all;
            is_first_frame = false;
            
            time_array = [time_array; 0];
            msd_array = [msd_array; 0];
            fprintf('Initialized Frame 1 (Timestep %d)\n', timestep);
            continue;
        end
        
        % Calculate raw step for ALL atoms
        dr = coords_all - r_prev;
        
        % Minimum Image Convention
        dr_corrected = dr - L * round(dr / L);
        
        % Accumulate the continuous unwrapped trajectory
        r_unwrap = r_unwrap + dr_corrected;
        
        % Calculate squared displacement from the origin (r0)
        disp_sq = sum((r_unwrap - r0).^2, 2);
        
        % Average over ALL atoms in the system
        current_msd = mean(disp_sq);
        
        current_time_ps = (timestep * dt_fs) / 1000;
        
        time_array = [time_array; current_time_ps];
        msd_array = [msd_array; current_msd];
        
        r_prev = coords_all;
        
        fprintf('Processed Timestep %d | Total MSD: %.4f A^2\n', timestep, current_msd);
    end
end
fclose(fid);

% =========================================================================
% Plotting and Total Diffusion Extraction
% =========================================================================
fig = figure('Visible','off');
plot(time_array, msd_array, 'k-', 'LineWidth', 2);
xlabel('Time (ps)', 'FontSize', 12);
ylabel('Total MSD (\AA^2)', 'FontSize', 12);
title('Total System Mean Square Displacement', 'FontSize', 14);
grid on;

filename = "msd" + ".png";
saveas(fig, filename);
close(fig);

% Extract D_total from the linear regime
fit_idx = floor(length(time_array)*0.2):length(time_array); 
p = polyfit(time_array(fit_idx), msd_array(fit_idx), 1);
slope = p(1); 

% D = slope / 6 (Convert A^2/ps to cm^2/s)
D_total_cm2_s = (slope / 6) * 1e-4;

fprintf('\n=======================================\n');
fprintf('TOTAL SYSTEM DIFFUSION ANALYSIS:\n');
fprintf('Slope of linear regime : %.4e A^2/ps\n', slope);
fprintf('Total Diffusion (D_tot): %.4e cm^2/s\n', D_total_cm2_s);
fprintf('=======================================\n');




% =========================================================================
% Viscosity Calculation via Stokes-Einstein Equation
% =========================================================================


% 2. Thermodynamic Constants
T_kelvin = 1733;                       % Simulation Temperature (K)
k_B = 1.380649e-23;                      % Boltzmann Constant (J/K or kg*m^2/(s^2*K))

% 3. Effective Ionic Radii (from paper, converted from Angstroms to meters)
% 1 Angstrom = 1e-10 meters
r_O_m  = 1.40 * 1e-10;  % O^2- radius (often used for bulk network flow)

% 4. Convert Diffusion Coefficients from cm^2/s to m^2/s
% 1 cm^2 = 1e-4 m^2
D_total_m2_s = D_total_cm2_s * 1e-4;

% =========================================================================
% 5. Execute Stokes-Einstein Equation: eta = (k_B * T) / (6 * pi * r * D)
% =========================================================================

% A. Total System Viscosity (The physically meaningful bulk property)
% Using the Oxygen radius as it dictates the primary structural network
eta_total_Pas = (k_B * T_kelvin) / (6 * pi * r_O_m * D_total_m2_s);


% Convert to Poise (optional, heavily used in older metallurgy texts)
% 1 Pa.s = 10 Poise (P)
eta_total_poise = eta_total_Pas * 10;

% =========================================================================
% 6. Output Results
% =========================================================================
fprintf('\n======================================================\n');
fprintf('VISCOSITY RESULTS (Stokes-Einstein at %d K):\n', T_kelvin);
fprintf('======================================================\n');
fprintf('Total System Viscosity (Bulk Flow):\n');
fprintf('  eta = %.4f Pa.s\n', eta_total_Pas);
fprintf('  eta = %.4f Poise\n\n', eta_total_poise);
