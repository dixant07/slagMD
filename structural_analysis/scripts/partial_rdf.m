clear all;
clc;

filename = 'dump.random_full.lammpstrj'; 
L = 17.750683960000;                     
rmax = 8.87534198;
dr = 0.1;

trj = readlines(filename);
line = 1;
avg_Na=0;
avg_Nb=0;
% 1:Ca, 2:Si, 3:Mg, 4:Al, 5:O
typeA = 1; % central atom
typeB = 2; % target atom

nbins = floor(rmax / dr);
g_hist = zeros(nbins, 1);
frames = 0;
r_axis = (0.5:1:nbins)' * dr;

while line<length(trj)-1
    if trj(line) == "ITEM: TIMESTEP"
        timestep = str2double(trj(line + 1));
        line = line + 1;

    elseif trj(line) == "ITEM: NUMBER OF ATOMS"
        numAtoms = str2double(trj(line + 1));
        line = line + 1;
    
    elseif trj(line) == "ITEM: ATOMS id type q x y z"
        atomData=str2double(split(trj(line+1:line+numAtoms), ' '));
        type = atomData(:,2); % atom type
        coord = atomData(:,4:6); % atom coordinates
        coordA = coord(type==typeA, :);
        coordB = coord(type==typeB, :);
        
        Na=length(coordA);
        Nb=length(coordB);

        avg_Na = avg_Na + Na;
        avg_Nb = avg_Nb + Nb;
        
        for i = 1:Na
            if typeA == typeB
                % for avoid double counting
                delta = coordB((i+1):end, :) - coordA(i, :);
                add_val = 2;
            else
                delta = coordB - coordA(i, :);
                add_val = 1;
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
        line = line+numAtoms;
        frames = frames+1;
    else
        line=line+1;
    end
end


avg_Na = avg_Na / frames;
avg_Nb = avg_Nb / frames;

rho_B = avg_Nb / (L^3);


g_AB = zeros(nbins, 1);
for k = 1:nbins
    r_inner = (k - 1) * dr;
    r_outer = k * dr;
    
    V_shell = (4/3) * pi * (r_outer^3 - r_inner^3);
    n_ideal_B = rho_B * V_shell;
    
    g_AB(k) = g_hist(k) / (frames * avg_Na * n_ideal_B);
end


fig = figure(1);

plot(r_axis, g_AB, 'r-', 'LineWidth', 3.5);
xlabel('Distance r (Angstroms)', 'FontSize', 16);
ylabel('g(r)', 'FontSize', 16);
title('Al-O RDF', 'FontSize', 18);
grid on;
xlim([0 rmax]);

[gmax, idx_max] = max(g_AB);
r_at_max = r_axis(idx_max);

for i=idx_max:length(g_AB)
    if g_AB(i)<=g_AB(i+1)
        idx_min=i;
        break;
    else
        continue;
    end
end

r_first_min = r_axis(idx_min);

fprintf("Rcut for %d: %.2f\n", typeA, r_first_min);