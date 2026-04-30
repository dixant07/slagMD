clear all;
clc;

filename = 'dump.random_full.lammpstrj';

L = 44.0;
dt = 1.0; % fs
stride = 1000;

trj = readlines(filename);

time_arr = [];
msd_arr = [];

% target_type = 2; %for atom wise diffusivity uncomment this & set target
first = 1;
line = 1;

while line < length(trj)-1
    
    if trj(line) == "ITEM: TIMESTEP"
        
        timestep = str2double(trj(line+1));
        line = line + 1;
        
    elseif trj(line) == "ITEM: NUMBER OF ATOMS"
        
        N = str2double(trj(line+1));
        line = line + 1;
        
    elseif trj(line) == "ITEM: ATOMS id type q x y z"
        atomData = str2double(split(trj(line+1:line+N),' '));
        
        type = atomData(:,2);
        ids = atomData(:,1);
        typeIds = ids; % for atom wise diffusivity comment this
        % typeIds = ids(type==target_type,:); %for atom wise diffusivity uncomment
        coords = atomData(:,4:6);
        
        % sorting
        [~,ind] = sort(typeIds);
        coords = coords(ind,:);
        
        if first == 1
            r0 = coords;
            r_prev = coords;
            r_unwrap = coords;
            
            time_arr = [time_arr; 0];
            msd_arr = [msd_arr; 0];
            
            first = 0;
        else
            dr = coords - r_prev;
            dr = dr - L * round(dr./L);
            r_unwrap = r_unwrap + dr;
            diff = r_unwrap - r0;
            sq = sum(diff.^2,2);
            msd_now = mean(sq);
            t_ps = timestep * dt / 1000;
            time_arr = [time_arr; t_ps];
            msd_arr = [msd_arr; msd_now];
            r_prev = coords;
        end
        line = line + N;
    else
        line = line + 1;
    end
end


% plot
figure;
plot(time_arr, msd_arr, 'k-', 'LineWidth', 2);
xlabel('time (ps)');
ylabel('msd');
grid on;

% linear fit
start_i = floor(length(time_arr)*0.2);
p = polyfit(time_arr(start_i:end), msd_arr(start_i:end), 1);
slope = p(1);

D = slope / 6;
D = D * 1e-4; % cm2/s

fprintf('Slope = %f\n', slope);
fprintf('D = %e cm2/s\n', D);


% ================= viscosity =====================
T = 1823; %temp
kB = 1.380649e-23;
rO = 1.40e-10; %ionic radius
D_m = D * 1e-4;
eta = (kB * T) / (6 * pi * rO * D_m);
eta_poise = eta * 10;

fprintf('Viscosity Pa.s = %f\n', eta);
fprintf('Viscosity Poise = %f\n', eta_poise);