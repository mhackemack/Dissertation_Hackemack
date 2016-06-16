function phi=compu_phi(psi)
% compute scalar flux and current

global snq npar

% initialize
% % phi=zeros(npar.ndof,2);
phi=zeros(npar.ndof,1);

% compute moments
for idir=1:snq.sn
    i1=(idir-1)*npar.ndof + 1;
    i2=(idir  )*npar.ndof   ;
    phi = phi + snq.w(idir)* psi(i1:i2);
    %     phi = phi + SNQ.w(idir)* kron( psi(i1:i2),  [ 1 SNQ.mu(idir)] );
end
