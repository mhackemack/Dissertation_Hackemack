function build_matrices()
global mats npar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute elementary matrices
% mt  = mass matrix for total xs
% mt  = mass matrix for absorption xs
% ms  = mass matrix for scattering xs
% g   = gradient matrix
% d   = diffusion stiffness matrix
% e   = edge matrix e{1} for mu>0, e{2} for mu<0
% egr = edge gradient matrix
[mt, ma, ms, g, k, e, egr, egrC, qva] = compute_elem1();
mats.elem.mt    = mt;
mats.elem.ma    = ma;
mats.elem.ms    = ms;
mats.elem.g     = g;
mats.elem.k     = k;
mats.elem.e     = e;
mats.elem.egr   = egr;
mats.elem.egrC  = egrC;
mats.global.qva = qva;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute global matrices
% T = transport matrix acting on angular fluxes
% S = scattering operator acting on angular fluxes
[T,S_psi_psi,S_phi_psi,S_phi_phi]=compute_T1();
mats.global.T = T;
mats.global.S_psi_psi = S_psi_psi;
mats.global.S_phi_psi = S_phi_psi;
mats.global.S_phi_phi = S_phi_phi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute diffusion matrix
% A = diffusion system matrix
if npar.perform_DSA
    if strcmpi(npar.DSA_type, 'mip')
        mats.global.A = DSA_MIP();
    elseif strcmpi(npar.DSA_type, 'm4s')
        mats.global.A = DSA_M4S();
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%