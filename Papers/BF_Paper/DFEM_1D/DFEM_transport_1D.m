function DFEM_transport_1D

% clear all;
close all; % closes all plotting figures
clc;       % clear console screen

global npar dat snq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=16;
if mod(sn,2)~=0
    error('sn must be even')
end
% load the angular quadrature
loadquadrature(sn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
dataID=8;
% load data
porder=2; % select spatial approx order (1 or 2)
[qsa]=loadmydata(dataID,porder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npar.lump = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T= transport matrix acting on angular fluxes
% S= scattering operator acting on angular fluxes
% F= transport matrix acting on scalar fluxes
% D= discrete to moment matrix
% Sigma= scattering matrix acting on moments
% M= moment to discrete matri
% note: M.Sigma.D = S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build new scheme matrices
[T,S,qva]=build_matrices();

% direct solve for the angular flux
psi = (T-S) \ (qva+qsa) ;
% postprocess to obtain the scalar flux
phi=compu_phi(psi);

figID=10;
myplot2(figID,phi,npar.porder,npar.dx,'b-');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('position','FontSize',12);
ylabel('Scalar flux','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Done !!!!');
