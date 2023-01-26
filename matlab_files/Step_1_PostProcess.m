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

% Plot data:
% =========================================================================
ii = 1;

figure('color','w')
plot(x_p{ii},'k.')

figure('color','w')
plot(v_p{ii}(:,1),'k.')