function [L,S]=compute_T1(mt,ms,gg,e)

global snq npar

% compute the sum of the SN weights
sw = snq.sw;

% dimension of the linear system
nel = npar.nel;
ndof= npar.ndof;
porder=npar.porder;
n = npar.ndof * snq.sn;
L = spalloc(n,n,snq.sn*nel*(porder+1)^2);

% build the action of all transport sweeps on the scattering term
ndir = snq.sn;
for idir=1:ndir
    m=snq.mu(idir); % shortcut
    % choose the edge contribution based on the sweeping order (L to R or R to L)
    if(m>0)
        edg=e{1};
    else
        edg=e{2};
    end

    i1 = (idir-1)*ndof+1;
    i2 = idir*ndof;

    % transport matrix for this direction
    % = mu*gradient + total mass + edge
    % note that the jacobian dx/2 is already included in the mass matrices
    % and that it does not appear in the gradient and edge matrices for 1D
    Ld = ( m*gg + mt ) + m*edg;

    % one way is to build the entire SN matrix L
    % and the scattering matrix S
    L(i1:i2,i1:i2) = Ld;
    S(i1:i2,:) = kron(snq.w/sw , ms) ;
%     if(length(ms)==2)
%         S(i1:i2,:) = S(i1:i2,:) + kron(...
%             snq.mu .* snq.w, snq.mu(idir)*ms{2} ) /sw;
%     end

end


% W  = kron(sparse(diag(snq.w)),speye(n));            % angular weighting
% SS = kron( [ ones(snq.sn,1)/sw , snq.mu'], speye(n) ); % spherical harmonics
% SS=SS/sw;

