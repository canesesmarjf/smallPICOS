% In this script, we read in the magnetic field profile from Step 1 and we
% create the following plasma profiles and save then into a HDF5 file:
% 1 - Electron density nm(x)
% 2 - Parallel ion temperature, Tixm
% 3 - Perp ion temperature, Tiym
% 4 - Parallel ion bulk flow, Ui_x

clear all
close all
clc

saveFig  = 1;
saveData = 1;

% Import magnetic field data:
% =========================================================================
% Subscript "m" stands for "mesh" defined quantity
inputProfile.xm = h5read('fields_IC.h5','/x');
inputProfile.Bm.value = h5read('fields_IC.h5','/Bx');

% Create electron density profile:
% =========================================================================
nm_min = 1e19;
nm_max = 5e19;
xm_min = inputProfile.xm(1);
xm_max = inputProfile.xm(end);

nm_sigma = (xm_max-xm_min)/10;
inputProfile.nm.value = nm_min + (nm_max-nm_min)*gaussian(inputProfile.xm,0,nm_sigma);

% Create electron temperature profile: 
% =========================================================================
Tem_min = 0.5e3;
Tem_max = 1.5e3;
Tem_sigma = (xm_max-xm_min)/10;
inputProfile.Tem.value = Tem_min + (Tem_max-Tem_min)*gaussian(inputProfile.xm,-0.7,Tem_sigma) + (Tem_max-Tem_min)*gaussian(inputProfile.xm,+0.7,Tem_sigma);

% Create parallel ion temperature profile: 
% =========================================================================
Tixm_min = 0.5e3;
Tixm_max = 1e3;
Tixm_sigma = (xm_max-xm_min)/15;
inputProfile.Tixm.value = Tixm_min + (Tixm_max-Tixm_min)*gaussian(inputProfile.xm,0,Tixm_sigma);

% Create perpendicular ion temperature profile: 
% =========================================================================
Tiym_min = 0.5e3;
Tiym_max = 1.5e3;
Tiym_sigma = (xm_max-xm_min)/10;
inputProfile.Tiym.value = Tiym_min + (Tiym_max-Tiym_min)*gaussian(inputProfile.xm,-1,Tixm_sigma) + (Tiym_max-Tiym_min)*gaussian(inputProfile.xm,+1,Tixm_sigma);

% Create parallel flow profile:
% =========================================================================
Uxm_max = C_s(1e3,2);
inputProfile.Uxm.value = Uxm_max*tanh(inputProfile.xm);

% Plot results:
% ========================================================================= 
% Electron profiles:
figure('color','w');
subplot(2,1,1)
hold on
box on
hTe(1) = plot(inputProfile.xm,inputProfile.nm.value,'k');
legendText{1} = ['$n_{e}$'];

set(gca,'fontSize',14,'fontName','Times')
set(hTe,'lineWidth',2)
hL = legend(hTe, legendText);
hL.Interpreter = 'latex';
hL.FontSize = 18;

ylim([0,1.2]*max(inputProfile.nm.value))
xlim([xm_min,xm_max])
xlabel('x [m]','Interpreter','latex','FontSize',16)
ylabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',16)
title('Electron density','Interpreter','latex','FontSize',16)

subplot(2,1,2)
hold on
box on
hTe(1) = plot(inputProfile.xm,inputProfile.Tem.value,'k');
legendText{1} = ['$T_{e}$'];

set(gca,'fontSize',14,'fontName','Times')
set(hTe,'lineWidth',2)
hL = legend(hTe, legendText);
hL.Interpreter = 'latex';
hL.FontSize = 18;

ylim([0,1.2]*max(inputProfile.Tem.value))
xlim([xm_min,xm_max])
xlabel('x [m]','Interpreter','latex','FontSize',16)
ylabel('$T_e$ [eV]','Interpreter','latex','FontSize',16)
title('Electron temperature','Interpreter','latex','FontSize',16)

% Save figure:
if saveFig
    folderName = '';
    figureName = 'Step_2_ElectronProfiles';
    
    % PDF figure:
    exportgraphics(gcf,[folderName,figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600) 
end

% Ion profiles:
figure('color','w');
subplot(2,1,1)
hold on
box on
hT(1) = plot(inputProfile.xm,inputProfile.Tixm.value,'k');
legendText{1} = ['$T_{i\parallel}$'];
hT(2) = plot(inputProfile.xm,inputProfile.Tiym.value,'r');
legendText{2} = ['$T_{i\perp}$'];

set(gca,'fontSize',14,'fontName','Times')
set(hT,'lineWidth',2)
hL = legend(hT, legendText);
hL.Interpreter = 'latex';
hL.FontSize = 18;

ylim([0,1.2]*max(inputProfile.Tiym.value))
xlim([xm_min,xm_max])
xlabel('x [m]','Interpreter','latex','FontSize',16)
ylabel('$T_i$ [m]','Interpreter','latex','FontSize',16)
title('Ion temperature','Interpreter','latex','FontSize',16)

subplot(2,1,2)
hold on
box on
clear hT
clear legendText
hT(1) = plot(inputProfile.xm,inputProfile.Uxm.value,'k');
legendText{1} = ['$U_{\parallel}$'];

set(gca,'fontSize',14,'fontName','Times')
set(hT,'lineWidth',2)
hL = legend(hT, legendText);
hL.Interpreter = 'latex';
hL.FontSize = 18;

ylim(1.2*[-1,+1]*max(inputProfile.Uxm.value))
xlim([xm_min,xm_max])
xlabel('x [m]','Interpreter','latex','FontSize',16)
ylabel('$U_x$ [m/s]','Interpreter','latex','FontSize',16)
title('Plasma bulk flow velocity','Interpreter','latex','FontSize',16)

% Save figure:
if saveFig
    folderName = '';
    figureName = 'Step_2_IonProfiles';
    
    % PDF figure:
    exportgraphics(gcf,[folderName,figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600) 
end

% Create computational particles density:
% =========================================================================
nm_min = 0.3;
nm_max = 1;

nm_sigma = (xm_max-xm_min)/10;
inputProfile.nm_cp.value = nm_min + (nm_max-nm_min)*gaussian(inputProfile.xm,0,nm_sigma);

xx = inputProfile.xm;
yy = tanh((xx + 0.5)*7) - tanh((xx - 0.5)*7);
y = nm_min + (nm_max-nm_min)*yy/max(yy);

% Convert to PDF:
nm_cp_pdf = y/trapz(xx,y);

% Plot results:
figure('Color','w')
box on
hold on
plot(inputProfile.xm,nm_cp_pdf,'r','LineWidth',2)

set(gca,'fontSize',14,'fontName','Times')
set(hT,'lineWidth',2)
hL = legend(hT, legendText);
hL.Interpreter = 'latex';
hL.FontSize = 18;

ylim([0,1.2]*max(nm_cp_pdf))
xlim([xm_min,xm_max])
xlabel('x [m]','Interpreter','latex','FontSize',18)
ylabel('$n_m^{cp}$','Interpreter','latex','FontSize',18)
title('Computational particle density','Interpreter','latex','FontSize',16)

% Save data HDF5:
% =========================================================================
if saveData
    disp('Saving data ...')

    folderName = '';
    h5Name = 'electrons_IC.h5';
    fileName = [folderName,h5Name];
  
    try
        h5create(fileName,'/x'   ,size(inputProfile.xm));
        h5create(fileName,'/n'   ,size(inputProfile.xm));
        h5create(fileName,'/T' ,size(inputProfile.xm));     
    end
        h5write(fileName,'/x'   ,(inputProfile.xm));
        h5write(fileName,'/n'   ,(inputProfile.nm.value));
        h5write(fileName,'/T' ,(inputProfile.Tem.value));
    disp('Data saved!')
end

if saveData
    disp('Saving data ...')

    folderName = '';
    h5Name = 'ions_IC.h5';
    fileName = [folderName,h5Name];
  
    for ss = 1:3
        group = ['/ions_',num2str(ss)];
        try
            h5create(fileName,[group,'/x'],size(inputProfile.xm));           
            h5create(fileName,[group,'/Tpar'],size(inputProfile.xm));
            h5create(fileName,[group,'/Tper'],size(inputProfile.xm));
            h5create(fileName,[group,'/upar'],size(inputProfile.xm));     
        end
        h5write(fileName,[group,'/x'],(inputProfile.xm));
        h5write(fileName,[group,'/Tpar'],(inputProfile.Tixm.value));
        h5write(fileName,[group,'/Tper'],(inputProfile.Tiym.value));
        h5write(fileName,[group,'/upar'],(inputProfile.Uxm.value));
    end
    disp('Data saved!')
end

% Computational particle density:
if saveData
    disp('Saving data ...')

    folderName = '';
    h5Name = 'computationalParticles_IC.h5';
    fileName = [folderName,h5Name];
  
    try
        h5create(fileName,'/x'        ,size(inputProfile.xm));
        h5create(fileName,'/n_pdf',size(inputProfile.xm));
    end
        h5write(fileName,'/x'        ,(inputProfile.xm));
        h5write(fileName,'/n_pdf',(nm_cp_pdf));
    disp('Data saved!')
end

% Copy data to input_files directory: 
if saveData
! cp *.h5 ../input_files       
end
