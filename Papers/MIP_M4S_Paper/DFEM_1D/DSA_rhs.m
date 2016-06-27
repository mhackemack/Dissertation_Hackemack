function rhs = DSA_rhs(dphi)
% call on global structures
global mats
% build rhs vector
rhs = mats.global.S_phi_phi*dphi;