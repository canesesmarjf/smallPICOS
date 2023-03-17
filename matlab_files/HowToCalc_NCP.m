% understanding how to produce the number of computational particles per
% cell:
clear all
close all
clc

Nx = 50;
Nk = 8;
Lmin = -1;
Lmax = +1;
dx = (Lmax - Lmin)/Nx;
ii = 1:Nx;
xm = Lmin + (ii - 0.5)*dx;
pdf = @(y,dx) y/sum(y*dx);
ncp_profile =  gaussian(xm,0,0.6);

% Initial desired value:
N_CP = 15555;

% Initial estimate of ncp_m:
ncp_m = N_CP*pdf(ncp_profile,dx);

% Make sure that ncp_m is exactly divisable by dx and also multiple of Nk:
ncp_m = round(ncp_m*dx/Nk)*Nk/dx;

% Calculate new N_CP:
N_CP = sum(ncp_m)*dx;

disp(['N_CP    = ',num2str(N_CP)])
disp(['N_CP/Nk = ',num2str(N_CP/Nk)])

ncp_per_cell = ncp_m*dx;
ncp_per_cell(1:14)

figure; 
plot(xm,ncp_per_cell,'k.')
ylim([0,max(ncp_per_cell)*1.2])
