% Create magnetic field profiles:

clear all
close all

% Flags:
saveFig  = 1;
saveData = 1;

% Include magnetic field functions:
% =========================================================================
pathLibrary = "/Users/juanfcanesesmarin/Documents/MATLAB/MyRepos/MagneticFieldCode_Axisymmetric/";
addpath(genpath(pathLibrary));

% Dimensions of grid:
% =========================================================================
Nx = 2^8; % Axial coordinate
Nr = 1; % Radial coordinate

% Read coil geometry spreadsheet:
% =========================================================================
coilGeometry = readtable('Step_1_CoilGeometry.xlsx','Sheet','conf_1');

% Define coil currents per power supply:
% =========================================================================
cc = 1;
scenarioType = 1;
switch scenarioType
    case 1
        coilCurrents{cc}.Mirror      = 40*1e3; % Mirror coils
        coilCurrents{cc}.CentralCell = 4*1e3; % Central cell coils
    case 2
        coilCurrents{cc}.Mirror      = 50*1e3; % Mirror coils
        coilCurrents{cc}.CentralCell = 4*1e3; % Central cell coils
    case 3
        coilCurrents{cc}.Mirror      = 60*1e3; % Mirror coils
        coilCurrents{cc}.CentralCell = 4*1e3; % Central cell coils        
end

% Create "coils" structure:
% =========================================================================
% Based on "coilGeometry" and "coilCurrents", we produce an object called
% "coils" for each physical coil. The "coils" structure contains the
% position of each current filament loop, power supply current and current
% per loop.

disp('Creating "coils" structure...')
for cc = 1:numel(coilCurrents)
    [coils{cc}] = CreateCoilStructure(coilGeometry,coilCurrents{cc});
end
disp(['Complete!'])

% Plot Magnetic coils and current filaments:
% =========================================================================
if 1
    cc = 1;
    figure('color','w','Tag','coilSetup')
    hold on
    for ii = 1:numel(coils{cc})
        plot(coils{cc}{ii}.zfil,+coils{cc}{ii}.rfil,'r.');
        plot(coils{cc}{ii}.zfil,-coils{cc}{ii}.rfil,'r.');
    end
    % Formatting:
    set(gca,'FontName','times','FontSize',11)
    xlabel('x [m]','Interpreter','latex','FontSize',13)
    ylabel('r [m]','Interpreter','latex','FontSize',13)
    box on
    axis image
    xlim([coils{cc}{1}.z,coils{cc}{end}.z]*1.5)
    grid on

    if saveFig
        folderName = '';
        figureName = 'Step_1_CoilSetup';
        
        % PDF figure:
        exportgraphics(gcf,[folderName,figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 
    
        % TIFF figure:
        exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600) 
    end

end

% Define computational grid:
% =========================================================================
r1D = 1e-3;
x1D = linspace(-2  ,2  ,Nx)';

% Calculate magnetic field:
% =========================================================================
disp('Calculating magnetic field...')
t1 = tic;
for cc = 1:numel(coilCurrents)
    [Br{cc},Bx{cc},~,~,~,~] = CalculateMagField(coils{cc},x1D,r1D,'grid');
end
disp(['Complete! Elapsed time: ',num2str(toc(t1)),' s'])

% Magnetic field magnitude:
% =========================================================================
for cc = 1:numel(coilCurrents)
    B1D{cc} = sqrt(Br{cc}.*Br{cc}  + Bx{cc}.*Bx{cc});
end

% Plot magnetic field on-axis:
% =========================================================================
figure('color','w')
box on
hold on
lineColor = {'k','bl','r'};
for cc = 1:numel(coilCurrents)
    plot(x1D,B1D{cc}(:,1),lineColor{cc},'lineWidth',2)
    bMax(cc) = max(B1D{cc}(:,1));
end
ylim([0,1.2]*max(bMax));
set(gca,'FontName','times','FontSize',11)
xlabel('x [m]','Interpreter','latex','FontSize',13)
ylabel('B [T]','Interpreter','latex','FontSize',13)

if saveFig
    folderName = '';
    figureName = 'Step_1_MagneticField';
    
    % PDF figure:
    exportgraphics(gcf,[folderName,figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600) 
end


% Save data HDF5:
% =========================================================================
if saveData
    disp('Saving data ...')

    fileName = 'fields_IC.h5';
  
    try
        h5create(fileName,'/x',size(x1D));
        h5create(fileName,'/Bx',size(x1D));
        h5create(fileName,'/Ex',size(x1D));
    end

    h5write(fileName,'/x',x1D);
    h5write(fileName,'/Bx',B1D{cc}(:,1));
    h5write(fileName,'/Ex',zeros(size(x1D)));

    disp('Data saved!')
end