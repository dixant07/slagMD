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

        fgetl(fid); % ITEM: NUMBER OF ATOMS
        N = str2double(fgetl(fid));

        fgetl(fid); % ITEM: BOX BOUNDS
        box = zeros(3,2);
        for i = 1:3
            box(i,:) = sscanf(fgetl(fid), '%f %f');
        end

        fgetl(fid); % ITEM: ATOMS ...

        data = fscanf(fid, '%f %f %f %f %f %f', [6 N])';
        % columns: id type q x y z

        ids = data(:,1);
        pos = data(:,4:6);

        % sort by atom id to ensure consistency
        [ids, idx] = sort(ids);
        pos = pos(idx,:);

        if frame == 1
            r0 = pos; % reference positions
        end

        dr = pos - r0;
        dr2 = sum(dr.^2, 2);

        MSD(frame) = mean(dr2);
    end
end

fclose(fid);

% normalize time (optional)
t0 = timesteps(1);
time = timesteps - t0;

% plot
figure;
plot(time, MSD, 'LineWidth', 2);
xlabel('Time');
ylabel('MSD');
title('Mean Square Displacement');
grid on;
