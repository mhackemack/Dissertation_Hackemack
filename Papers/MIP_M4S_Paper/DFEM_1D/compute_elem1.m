function [mtot,mabs,msca,gr,e,qva]=compute_elem1()

global dat npar snq

% shortcuts
nel = npar.nel;
ndof=npar.ndof;
porder=npar.porder;
dx=npar.dx;

% enforce standard upwind scheme
fu=1; fd=0;
[xq,wq] = GLNodeWt(porder+1);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
g=m;
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

% compute local matrices + load vector
for i=1:porder+1
    for j=1:porder+1
        m(i,j)=  dot(wq.*b(:,i)    , b(:,j));
        g(i,j)= -dot(wq.*dbdx(:,i) , b(:,j));
        k(i,j)=  dot(wq.*dbdx(:,i) , dbdx(:,j));
    end
    f(i)= dot(wq, b(:,i));
end
if npar.lump
    m=diag(sum(m));
end
% jacobian of the affine transformation
jac  = dx(1:nel)/2;

% loop over elements
tot=zeros(nel,1); sca=tot; sabs=tot; qext=tot; DC=tot;
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    tot(iel)  = dat.sigt(my_zone);
    sabs(iel) = dat.siga(my_zone);
    sca(iel)  = dat.sigs(my_zone);
    DC(iel)   = dat.D(my_zone);
    qext(iel) = dat.qv(my_zone)/snq.sw; % normalize to sum weights, isotropic qext
end

% global mass matrices
mtot =  kron( sparse( diag(tot.*jac) ), sparse(m) );
mabs =  kron( sparse( diag(tot.*jac) ), sparse(m) );
msca =  kron( sparse( diag(sca.*jac) ), sparse(m) );
% global gradient matrix
gr = kron( speye(nel), sparse(g) );
% global stiffness matrix

% global exeternal source vector
qva = kron(qext.*jac,f);
% make it an angular ext source
qva=kron(ones(snq.sn,1),qva);

% edge matrix
e{1} = spalloc(ndof,ndof,ndof);
e{2} = spalloc(ndof,ndof,ndof);

% build edge matrix
for iel=1:nel
    % starting dofs indices for iel
    istart = (iel-1)*(porder+1)+1;
    % starting dofs indices for iel-1
    iel0 = iel-1;
    istart0 = (iel0-1)*(porder+1)+1;
    % starting dofs indices for iel+1
    iel2 = iel+1;
    istart2 = (iel2-1)*(porder+1)+1;
    
    % edge matrix, mu>0
    % ~~~~~~~~~~~~~~~~~
    if(iel~=1),      e{1}(istart  ,istart0+porder) = -fu; end
    e{1}(istart+porder,istart+porder ) =  fu;
    
    % edge matrix, mu<0
    % ~~~~~~~~~~~~~~~~~
    e{2}(istart  ,istart   ) = -fu;
    if(iel~=nel), e{2}(istart+porder,istart2  ) =  fu; end
end


