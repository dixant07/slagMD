clear all;
clc;

filename = 'dump.slag.lammpstrj';

L = 44.0;

% types
type_Si = 2;
type_Al = 4;
type_O  = 5;

% cutoff
R_SiO = 2.35;
R_AlO = 2.55;

R2_SiO = R_SiO^2;
R2_AlO = R_AlO^2;


trj = readlines(filename);

line = 1;

Qn = zeros(1,6); % Q0 Q1 Q2 Q3 Q4 other
total_Si = 0;
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
        
        nSi = length(coord_Si);
        nAl = length(coord_Al);
        nO  = length(coord_O);
        
        % adjacency
        A_SiO = zeros(nO,nSi);
        A_AlO = zeros(nO,nAl);
        
        % build Si-O
        for j = 1:nSi
            d = coord_O - coord_Si(j,:);
            d = d - L * round(d./L);
            r2 = sum(d.^2,2);
            A_SiO(:,j) = r2 < R2_SiO;
        end
        
        % build Al-O
        for k = 1:nAl
            d = coord_O - coord_Al(k,:);
            d = d - L * round(d./L);
            r2 = sum(d.^2,2);
            A_AlO(:,k) = r2 < R2_AlO;
        end
        
        % oxygen coordination
        ZO = sum(A_SiO,2) + sum(A_AlO,2);
        
        BO = (ZO == 2);
        
        % Qn
        n_vec = A_SiO' * BO;
        
        for j = 1:nSi
            nval = n_vec(j);
            
            if nval >= 0 && nval <= 4
                Qn(nval+1) = Qn(nval+1) + 1;
            else
                Qn(6) = Qn(6) + 1;
            end
        end
        
        total_Si = total_Si + nSi;
        frames = frames + 1;
        
        fprintf('Frame %d done\n',frames);
        
        line = line + N;
        
    else
        line = line + 1;
    end
    
end


pct = (Qn / total_Si) * 100;

fprintf('\nResults:\n');
fprintf('Q0 = %.2f %%\n',pct(1));
fprintf('Q1 = %.2f %%\n',pct(2));
fprintf('Q2 = %.2f %%\n',pct(3));
fprintf('Q3 = %.2f %%\n',pct(4));
fprintf('Q4 = %.2f %%\n',pct(5));
fprintf('Other = %.2f %%\n',pct(6));