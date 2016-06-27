clear; close all; clc;
global npar dat snq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin user input
% ------------------------------------------------------------------------------
% select angular approx (must be an even number)
sn=16;
loadquadrature(sn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
dataID=9;
% FEM order
porder=1; % select spatial approx order (e.g. 1, 2, etc.)
% lumping bool
npar.lump = 0;
% load data
[qsa]=loadmydata(dataID, porder);
% iterative procedure properties
maxits = 1e4;
tol = 1e-8;
% DSA properties
perform_DSA = false;
DSA_type = 'M4S';
% End user input
% ------------------------------------------------------------------------------
% Build all matrices
% ------------------------------------------------------------------------------
[T, S_psi_psi, S_phi_psi, S_phi_phi, qva]=build_matrices();
% Compute reference solution
% ------------------------------------------------------------------------------
% psi = (T-S_psi_psi) \ (qva+qsa) ;
phi_ref=compu_phi((T-S_psi_psi) \ (qva+qsa));
figID=10; myplot2(figID,phi_ref,npar.porder,npar.dx,'b-');
% Plot reference solution
xlabel('position','FontSize',12);
ylabel('Scalar flux','FontSize',12);
% Compute iterative solution
% ------------------------------------------------------------------------------
