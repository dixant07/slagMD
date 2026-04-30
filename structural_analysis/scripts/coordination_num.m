clear all;
clc;

filename = 'dump.random_full.lammpstrj'; 
L = 17.750683960000;                     

trj = readlines(filename);

% 1:Ca, 2:Si, 3:Mg, 4:Al, 5:O
typeA = 1; %central atom
typeB = 2; %target atom

R_cut = 2.25; %from rdf

total_neighbors = 0;
total_A_atoms = 0;
frames = 0;
line = 1;

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
    
        for i=1:Na
            delta = coordB - coordA(i, :);
            delta = delta - L * round(delta / L);
            r = sqrt(sum(delta.^2, 2));
            neighbors = sum(r > 0 & r < R_cut);
            total_neighbors = total_neighbors + neighbors;
        end
        total_A_atoms = total_A_atoms + Na;
        line = line+numAtoms;
        frames = frames+1;
    else
        line=line+1;
    end

end

final_CN = total_neighbors / total_A_atoms;

fprintf('\n=======================================\n');
fprintf('FINAL COORDINATION NUMBER RESULTS:\n');
fprintf('Central Atom Type : %d\n', typeA);
fprintf('Neighbor Atom Type: %d\n', typeB);
fprintf('Cutoff Radius     : %.2f A\n', R_cut);
fprintf('Average CN        : %.4f\n', final_CN);
fprintf('=======================================\n');
