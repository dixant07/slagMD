clear; clc;

filename = 'dump.slag.lammpstrj';
fid = fopen(filename, 'r');

timesteps = [];
MSD = [];

frame = 0;

while ~feof(fid)

    line = fgetl(fid);

    if contains(line, 'ITEM: TIMESTEP')
        frame = frame + 1;

        t = str2double(fgetl(fid));
        timesteps(frame) = t;

        fgetl(fid); % NUMBER OF ATOMS
        N = str2double(fgetl(fid));

        fgetl(fid); % BOX BOUNDS
        box = zeros(3,2);
        for i = 1:3
            box(i,:) = sscanf(fgetl(fid), '%f %f');
        end

        L = box(:,2) - box(:,1); % box length [Lx Ly Lz]

        fgetl(fid); % ATOMS header

        data = fscanf(fid, '%f %f %f %f %f %f', [6 N])';
        ids = data(:,1);
        pos = data(:,4:6);

        % sort atoms by id
        [ids, idx] = sort(ids);
        pos = pos(idx,:);

        if frame == 1
            r0 = pos;
            r_prev = pos;

            unwrap = zeros(N,3); % cumulative correction
            r_unwrapped = pos;

        else
            dr = pos - r_prev;

            % periodic correction (minimum image)
            for d = 1:3
                over = dr(:,d) >  L(d)/2;
                under = dr(:,d) < -L(d)/2;

                unwrap(over,d)  = unwrap(over,d)  - L(d);
                unwrap(under,d) = unwrap(under,d) + L(d);
            end

            r_unwrapped = pos + unwrap;
            r_prev = pos;
        end

        % MSD using unwrapped coordinates
        dr_total = r_unwrapped - r0;
        MSD(frame) = mean(sum(dr_total.^2, 2));
    end
end

fclose(fid);

time = timesteps - timesteps(1);

figure;
plot(time, MSD, 'LineWidth', 2);
xlabel('Time');
ylabel('MSD');
title('Mean Square Displacement');
grid on;


% ===== SELECT LINEAR REGION =====
start_idx = round(0.4 * length(time));   % adjust based on your plot
end_idx   = round(0.9 * length(time));

t_fit = time(start_idx:end_idx);
msd_fit = MSD(start_idx:end_idx);

p = polyfit(t_fit, msd_fit, 1);
slope = p(1);   % d(MSD)/dt

D = slope / 6;   % Å^2/fs

% ===== INPUTS =====
T = 1823.15;                  % temperature (K)
d_A = 3.0;                 % effective diameter (Å)  <-- must choose carefully
c = 6;                     % 6 = stick, 4 = slip

% ===== CONSTANTS =====
kB = 1.380649e-23;         % J/K

% ===== UNIT CONVERSION =====
% Å^2/fs -> m^2/s
D = D * 1e-5;

% Å -> m
d = d_A * 1e-10;

% ===== VISCOSITY =====
eta = kB * T / (c * pi * D * d);

% ===== OUTPUT =====
fprintf('Diffusivity (SI): %.3e m^2/s\n', D);
fprintf('Viscosity: %.3e Pa.s\n', eta);