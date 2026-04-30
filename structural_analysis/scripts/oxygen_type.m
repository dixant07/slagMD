clear all;
clc;

filename = 'dump.slag.lammpstrj';

L = 44.0;

% types
type_Si = 2;
type_Al = 4;
type_O  = 5;

% cutoff (from RDF)
R_SiO = 2.35;
R_AlO = 2.55;

trj = readlines(filename);

line = 1;

sum_FO = 0;
sum_NBO = 0;
sum_BO = 0;
sum_TRI = 0;

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
        
        coord_Si = coord(type==type_Si,:);
        coord_Al = coord(type==type_Al,:);
        coord_O  = coord(type==type_O,:);
        
        nO = length(coord_O);
        
        c_FO = 0;
        c_NBO = 0;
        c_BO = 0;
        c_TRI = 0;
        
        for i = 1:nO
            
            % Si bonds
            d_Si = coord_Si - coord_O(i,:);
            d_Si = d_Si - L * round(d_Si./L);
            r_Si = sqrt(sum(d_Si.^2,2));
            b_Si = sum(r_Si < R_SiO);
            
            % Al bonds
            d_Al = coord_Al - coord_O(i,:);
            d_Al = d_Al - L * round(d_Al./L);
            r_Al = sqrt(sum(d_Al.^2,2));
            b_Al = sum(r_Al < R_AlO);
            
            Z = b_Si + b_Al;
            
            if Z == 0
                c_FO = c_FO + 1;
            elseif Z == 1
                c_NBO = c_NBO + 1;
            elseif Z == 2
                c_BO = c_BO + 1;
            else
                c_TRI = c_TRI + 1;
            end
            
        end
        
        pct_FO = (c_FO / nO) * 100;
        pct_NBO = (c_NBO / nO) * 100;
        pct_BO = (c_BO / nO) * 100;
        pct_TRI = (c_TRI / nO) * 100;
        
        sum_FO = sum_FO + pct_FO;
        sum_NBO = sum_NBO + pct_NBO;
        sum_BO = sum_BO + pct_BO;
        sum_TRI = sum_TRI + pct_TRI;
        
        frames = frames + 1;
        
        fprintf('Frame %d | BO %.2f | NBO %.2f | FO %.2f\n',frames,pct_BO,pct_NBO,pct_FO);
        
        line = line + N;
        
    else
        line = line + 1;
    end
    
end


final_FO = sum_FO / frames;
final_NBO = sum_NBO / frames;
final_BO = sum_BO / frames;
final_TRI = sum_TRI / frames;

fprintf('\nFinal Results:\n');
fprintf('FO = %.2f %%\n',final_FO);
fprintf('NBO = %.2f %%\n',final_NBO);
fprintf('BO = %.2f %%\n',final_BO);
fprintf('TRI = %.2f %%\n',final_TRI);