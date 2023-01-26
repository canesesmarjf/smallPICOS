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

% Calculate histograms for each species:
% =========================================================================
x_m = d.fields.x_m;
dx = x_m(2) - x_m(1);
x_edges = x_m(1:end-1) + dx/2;
[counts_0,~] = histcounts(d.ions_0.x_p,x_edges,"Normalization","count");
[counts_1,~] = histcounts(d.ions_1.x_p,x_edges,"Normalization","count");

% Plot data:
% =========================================================================
figure('color','w')
subplot(2,1,1)
box on
hold on
plot(d.fields.x_m,d.ions_0.ncp_m*dx,'k.-')
plot(d.fields.x_m,d.ions_1.ncp_m*dx,'r.-')

plot(x_m(2:end-1),counts_0,'k','LineWidth',2)
plot(x_m(2:end-1),counts_1,'r','LineWidth',2)

ylim([0,1.2]*max(d.ions_0.ncp_m*dx))

subplot(2,1,2)
box on
hold on
plot(d.fields.x_m,d.ions_0.n_m,'k.-')
plot(d.fields.x_m,d.ions_1.n_m,'r.-')
ylim([0,1.2]*max(d.ions_0.n_m))

figure;
plot(d.ions_0.x_p,d.ions_0.v_p(:,1),'k.')

if 0
    figure;
    plot3(d.ions_0.x_p,d.ions_0.v_p(:,1),d.ions_0.v_p(:,2),'k.')
end

figure;
plot(d.fields.x_m,d.fields.Bx_m,'k.-')