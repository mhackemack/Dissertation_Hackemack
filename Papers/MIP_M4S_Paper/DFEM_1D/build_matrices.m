function [T, S_psi_psi,S_phi_psi,S_phi_phi, qva]=build_matrices()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute elementary matrices
% mt = mass matrix for total xs
% ms = mass matrix for scattering xs
% g  = gradient matrix
% e =  edge matrix e{1} for mu>0, e{2} for mu<0
[mt , ma, ms ,g ,e,qva] = compute_elem1();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute global matrices
% T= transport matrix acting on angular fluxes
% S= scattering operator acting on angular fluxes
% F= transport matrix acting on scalar fluxes
% D= discrete to moment matrix
% Sigma= scattering matrix acting on moments
% M= moment to discrete matrix
% note: M.Sigma.D = S
% C= matrix acting on angular fluxes and producing the current

[T,S_psi_psi,S_phi_psi,S_phi_phi]=compute_T1(mt,ma,ms,g,e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

