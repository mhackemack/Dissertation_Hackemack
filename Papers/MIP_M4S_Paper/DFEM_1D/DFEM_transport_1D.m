clear; close all; clc;
global npar dat snq mats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin user input
% ------------------------------------------------------------------------------
% select angular approx (must be an even number)
sn=8;
loadquadrature(sn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
dataID=10;
% FEM order
porder=1; % select spatial approx order (e.g. 1, 2, etc.)
% lumping bool
npar.lump = 0;
% load data
loadmydata(dataID, porder);
% iterative procedure properties
npar.maxits = 1e6;
npar.tol = 1e-8;
% DSA properties
npar.perform_DSA = 1;
npar.DSA_type    = 'MIP';
% plotting options
plot_bool = false;
% End user input
% ------------------------------------------------------------------------------
% Build all matrices
% ------------------------------------------------------------------------------
build_matrices();
psi_tot_src = mats.global.qva+mats.global.qsa;
% Compute reference solution
% ------------------------------------------------------------------------------
phi_ref=compu_phi((mats.global.T-mats.global.S_psi_psi) \ (psi_tot_src));
% Plot reference solution
if plot_bool
    figID = subplot(2,1,1);
    pref = myplot2(figID,phi_ref,npar.porder,npar.dx,'b-');
    xlabel('position','FontSize',12);
    ylabel('Scalar flux','FontSize',12);
end
% Compute iterative solution
% ------------------------------------------------------------------------------
phi = zeros(npar.ndof,1); phi0 = phi;
count = 0; err0 = 0;
% Loop through iterations
for i=1:npar.maxits
    % update counter
    count = count + 1;
    % perform transport calculation
    psi = (mats.global.T)\(psi_tot_src + mats.global.S_phi_psi*phi);
    phi = compu_phi(psi);
    % perform dsa step
    if npar.perform_DSA
        dphi = (mats.global.A)\DSA_rhs(phi-phi0);
        phi = phi + dphi;
    end
    % check residual
    err = norm(phi-phi0)/norm(phi);
    if i==1
        msg = sprintf('Iteration %d - error = %e',i,err);
    else
        msg = sprintf('Iteration %d - error = %e, NSR = %e',i,err,err/err0);
    end
    disp(msg);
    if err < npar.tol
        break;
    else
        phi0 = phi;
        err0 = err;
    end
end
% Print difference between reference and iterative solutions
fprintf('\nAbsolute Difference = %e\n',norm(phi-phi_ref))
fprintf('Relative Difference = %e\n',norm(phi-phi_ref)/norm(phi))
if plot_bool
    % plot iterative solution
    pact = myplot2(figID,phi,npar.porder,npar.dx,'r--');
    legend([pref,pact],{'\phi_{ref}','\phi_h'},'Location','Best');
    % plot difference between iterative and reference solutions
    figID = subplot(2,1,2);
    myplot2(figID,abs(phi-phi_ref),npar.porder,npar.dx,'k-');
    xlabel('position','FontSize',12);
    ylabel('|\phi_h - \phi_{ref}|','FontSize',12);
end