clear all;
clc;

filename = 'dump.slag.lammpstrj';

L = 44.0;

% types
type_Al = 4;
type_O  = 5;

% cutoff
R_AlO = 2.55;
R2 = R_AlO^2; % square cutoff

trj = readlines(filename);

line = 1;

count_AlO4 = 0;
count_AlO5 = 0;
count_other = 0;

total_Al = 0;
frames = 0;

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
        coord = atomData(:,4:6);
        
        coord_Al = coord(type==type_Al,:);
        coord_O  = coord(type==type_O,:);
        
        nAl = length(coord_Al);
        
        c4 = 0;
        c5 = 0;
        cO = 0;
        
        for k = 1:nAl
            
            d = coord_O - coord_Al(k,:);
            
            % minimum image
            d = d - L * round(d./L);
            
            r2 = sum(d.^2,2);
            
            Z = sum(r2 < R2);
            
            if Z == 4
                c4 = c4 + 1;
            elseif Z == 5
                c5 = c5 + 1;
            else
                cO = cO + 1;
            end
            
        end
        
        count_AlO4 = count_AlO4 + c4;
        count_AlO5 = count_AlO5 + c5;
        count_other = count_other + cO;
        
        total_Al = total_Al + nAl;
        frames = frames + 1;
        
        fprintf('Frame %d | AlO4 %d | AlO5 %d\n',frames,c4,c5);
        
        line = line + N;
        
    else
        line = line + 1;
    end
    
end


pct_AlO4 = (count_AlO4 / total_Al) * 100;
pct_AlO5 = (count_AlO5 / total_Al) * 100;
pct_other = (count_other / total_Al) * 100;

fprintf('\nFinal:\n');
fprintf('AlO4 = %.2f %%\n',pct_AlO4);
fprintf('AlO5 = %.2f %%\n',pct_AlO5);
fprintf('Other = %.2f %%\n',pct_other);