% This script was made to read data produced by smallPICOS++:
clear all
close all

% 1 - Define file names and path:
% =========================================================================
fileName = 'file_1.h5';
dirName = '../output_files/HDF5/';
pathName = [dirName,fileName];

% 2- Get group names:
% =========================================================================
info = h5info(pathName);
for i = 1:numel(info.Groups)
    groupNames{i} = info.Groups(i).Name;
    for j = 1:numel(info.Groups(i).Datasets)
        datasetNames{i}{j} = info.Groups(i).Datasets(j).Name;
        info.Groups(i).Datasets(j).Name
    end
end

% 3 - Read data:
% =========================================================================
for i = 1:numel(info.Groups)
    group = groupNames{i}(2:end);
    disp(["Group = ",group]);
    for j = 1:numel(info.Groups(i).Datasets)
        var = info.Groups(i).Datasets(j).Name;
        values = h5read(pathName,[groupNames{i},'/',datasetNames{i}{j}]);
        
        d.(group).(var) = values;
        data{i}{j} = values;
    end
end 

% Assign data:
% =========================================================================
for i = 1:numel(info.Groups)
    a_p{i} = data{i}{1};
    v_p{i} = data{i}{2};
    x_p{i} = data{i}{3};
end 


figure('color','w')
subplot(2,1,1)
box on
hold on
x_m = d.fields.x_m;
dx = x_m(2) - x_m(1);
plot(d.fields.x_m,d.ions_0.ncp_m*dx,'k.-')
plot(d.fields.x_m,d.ions_1.ncp_m*dx,'r.-')
ylim([0,1.2]*max(d.ions_0.ncp_m*dx))

subplot(2,1,2)
box on
hold on
plot(d.fields.x_m,d.ions_0.n_m,'k.-')
plot(d.fields.x_m,d.ions_1.n_m,'r.-')
ylim([0,1.2]*max(d.ions_0.n_m))

figure;
plot(d.ions_0.x_p,d.ions_0.v_p(:,1),'k.')

figure;
plot3(d.ions_0.x_p,d.ions_0.v_p(:,1),d.ions_0.v_p(:,2),'k.')

figure;
plot(d.fields.x_m,d.fields.Bx_m,'k.-')

return

% Plot data:
% =========================================================================

figure('color','w')
hold on
for ii = 1:numel(groupNames)
    if ~strcmpi(groupNames{ii},'x_p')
        continue
    end
    hx(ii) = plot(x_p{ii});
    set(hx(ii),'lineStyle','none','marker','.');
end

figure('color','w')
hold on
for ii = 1:numel(groupNames)
    hv(ii) = plot(x_p{ii},v_p{ii}(:,1));
    set(hv(ii),'lineStyle','none','marker','.');
end