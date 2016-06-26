function F=hc_ss_x_CG_gap
% Solves the time-dependent heat conduction equation in 1-D x-geometry 
% using CFEM with T gap.
% The conductivities and the volumetric sources can be spatially dependent.

% clear the console screen
clc; clear all; close all;
% load the data structure with info pertaining to the physical problem
dat.k{1}=@k_Zr; % W/m-K
dat.k{2}=@k_fuel; % W/m-K
dat.k{3}=@k_clad; % W/m-K
dat.esrc{1}=@zero_function; % W/m3
dat.esrc{2}=@esrc; % W/m3
dat.esrc{3}=@zero_function; % W/m3

dat.hgap=15764; % W/m^2-C
dat.hcv=1612.414; % W/m^2-C
dat.width=[0.003175 0.0174115 0.0179195]; % m
bc.left.type=0; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
bc.rite.type=1;
bc.rite.C=50;
dat.bc=bc; clear bc;

gap_zone_ID=2;
nel_zone = [ 10 30 2];

% load the numerical parameters, npar, structure pertaining to numerics
% number of elements
if length(dat.width)~=length(nel_zone)
    error('not the same number of zones in dat.width and nel_zone');
end
if length(dat.k)~=length(nel_zone)
    error('not the same number of zones in dat.k and nel_zone');
end
if length(dat.esrc)~=length(nel_zone)
    error('not the same number of zones in dat.esrc and nel_zone');
end
npar.nel = sum(nel_zone);
npar.nelgap = sum(nel_zone(1:gap_zone_ID));

% domain
tmp_width = [ 0 dat.width];
x=[];
iel2zon=[];
for z=1:length(nel_zone)
    x_zone = linspace(tmp_width(z),tmp_width(z+1),nel_zone(z)+1);
    if isempty(x)
        x = x_zone;
        iel2zon=z*ones(nel_zone(z),1);
    else
        x = [x x_zone(2:end)];
        iel2zon =[ iel2zon; z*ones(nel_zone(z),1)];
    end
end
npar.x=x;
npar.iel2zon=iel2zon;
% polynomial degree
npar.porder=2;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+2;
% connectivity

gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
% adding the gap unknowns
gn(npar.nelgap+1:end,1:end)=gn(npar.nelgap+1:end,1:end)+1;
npar.gn=gn; clear gn;

% solve system
F = solve_fem3(dat,npar);
% plot
figure(1)

% create x values for plotting FEM solution
npar.xf=zeros(npar.nel*npar.porder+1,1);
for iel=1:npar.nel
    ibeg = (iel-1)*npar.porder+1;
    iend = (iel  )*npar.porder+1;
    npar.xf(ibeg:iend)=linspace(npar.x(iel),npar.x(iel+1),npar.porder+1) ;
end
npar.xfi=[npar.xf(1:npar.nelgap*npar.porder+1);npar.xf(npar.nelgap*npar.porder+1:end)] ;

% verification is always good
a=verif_hc_eq(dat);

% get coefficient for the analytical solution
k=dat.k; src=dat.esrc; L=dat.width;
x1=linspace(0,L(1));
x2=linspace(L(1),L(2));
x3=linspace(L(2),L(3));
y1=a(1)*x1+a(2);
y2=-src{2}(x2)/(2*k{2}(x2))*(x2.^2)+a(3)*x2+a(4);
y3=a(5)*x3+a(6);

plot(npar.xfi,F,'.-',x1,y1,'r-',x2,y2,'r-',x3,y3,'r-'); hold all;
title('1D steady-state heat conduction, with T gap, Cartesian coordinates')
xlabel('Width (m)')
ylabel('Temperature (C)')
legend('FEM','Analytical','Location','northoutside','Orientation','horizontal')

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=solve_fem3(dat, npar)

% initial guess
% (not needed if a direct method is used to solve the linear system)
u=ones(npar.ndofs,1);

% assemble the matrix and the rhs
[A,b]=assemble_system(npar,dat);

% solve
u=A\b;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs]=assemble_system(npar,dat)

% assemble the matrix, the rhs, apply BC

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
g1 = gn(npar.nelgap,end);
g2 = gn(npar.nelgap+1,1);
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+3)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*porder+2;
% allocate memory
A=spalloc(n,n,nnz);
rhs=zeros(n,1);

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
poly_max=2*porder;
[xq,wq] = GLNodeWt(porder+1);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      = int_domain (b rhs)

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    my_zone=npar.iel2zon(iel);
    d=dat.k{my_zone}(x);
    q=dat.esrc{my_zone}(x);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            % m(i,j)= dot(0.*wq.*b(:,i)    , b(:,j));
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    % A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + ...
    %    m*Jac + k/Jac;
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end

A(g1,g1)=A(g1,g1)+dat.hgap;
A(g1,g2)=A(g1,g2)-dat.hgap;
A(g2,g1)=A(g2,g1)-dat.hgap;
A(g2,g2)=A(g2,g2)+dat.hgap;

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+dat.hcv;
        rhs(1)=rhs(1)+dat.hcv*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+dat.hcv;
        rhs(n)=rhs(n)+dat.hcv*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
    rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shapefun,dhdx]=feshpln (xv,p)
% computes the finite element basis functions and their first derivatives
% for any order p
% here, we use Lagrangian FE functions

xd=linspace(-1,1,p+1);

shapefun=zeros(length(xv),p+1);
dhdx    =zeros(length(xv),p+1);

% shape function
for i=1:p+1
    num=1.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=num.*(xv-xd(j));
            den=den.*(xd(i)-xd(j));
        end
    end
    shapefun(:,i)=num./den;
end

% derivative of the shape function
for i=1:p+1
    sum=0.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=1;
            for k=1:p+1
                if((k~=i)&&(k~=j))
                    num=num.*(xv-xd(k));
                end
            end
            sum=sum+num;
            den=den.*(xd(i)-xd(j));
        end
    end
    dhdx(:,i)=sum./den;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = GLNodeWt(n)
% GLNodeWt  Nodes and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [x,w] = GLNodeWt(n)
%
% Input:     n = order of quadrature rule
%
% Output:    x = vector of nodes
%            w = vector of weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = 2*V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a=verif_hc_eq(dat)

k=dat.k; src=dat.esrc; hgap=dat.hgap; hcv=dat.hcv; L=dat.width;

% general form of the solution:
% Zone 1 : T1 = B1*x + E1
% dT1/dx= B1
% Zone 2 : T2 = -q/(2*k2)*(x.^2) + B2*x + E2
% dT2/dx= -q/k2*x + B2
% Zone 3 : T3 = B3*x + E3
% dT3/dx= B3

switch dat.bc.left.type
    case 0 % Neumann
        % k1*du1/dn=C on the left becomes: -k1*du1/dx=C
        % <==> -k1*B1=C <==> B1=-C/k1
        mat(1,1:6) =[1,0,0,0,0,0];
        b(1) = -dat.bc.left.C / k{1}(0);
    case 1 % Robin
        % u1+k1/hcv*du1/dn =C on the left becomes: u1-k1/hcv*du1/dx =C
        % <==> -k1/hcv*B1+E1=C
        mat(1,1:6) =[-k{1}(0)/hcv,1,0,0,0,0];
        b(1) = dat.bc.left.C;
    case 2 % Dirichlet
        % u1=C <==> E1=C
        mat(1,1:6) =[0,1,0,0,0,0];
        b(1) = dat.bc.left.C;
end
switch dat.bc.rite.type
    case 0 % Neumann
        % k3*du3/dn=C on the right becomes: k3*du3/dx=C
        % <==> k3*B3=C <==> B3=C/k3
        mat(6,1:6) =[0,0,0,0,1,0];
        b(6) = dat.bc.rite.C / k{3}(L(3));
    case 1 % Robin
        % u3+k3/hcv*du3/dn =C on the right becomes: u3+k3/hcv*du3/dx =C
        % <==> (L3+k3/hcv)B3+E3=C
        mat(6,1:6) =[0,0,0,0,L(3)+k{3}(L(3))/hcv,1];
        b(6) = dat.bc.rite.C;
    case 2 % Dirichlet
        % u3=C <==> B3*L3+E3=C
        mat(6,1:6) =[0,0,0,0,L(3),1];
        b(6) = dat.bc.rite.C;
end

% fixed conditions
% continuity of T and flux between zone 1 and zone 2 (interface L1)
% T1(L1)=T2(L1) <==> B1*L1+E1=(-q/2k2)*(L1^2)+B2*L1+E2
% <==> B1*L1+E1-B2*L1-E2=(-q/2k2)*(L1^2)
mat(2,1:6) =[L(1),1,-L(1),-1,0,0];
b(2) =-src{2}(L(1))/(2*k{2}(L(1)))*L(1)*L(1);
% phi1(L1)=phi2(L1) <==> k1*B1=k2*((-q/k2)*L1+B2)
% <==> (k1/k2)*B1-B2=(-q/k2)*L1
mat(3,1:6) =[k{1}(L(1))/k{2}(L(1)),0,-1,0,0,0];
b(3) =-src{2}(L(1))/k{2}(L(1))*L(1);

% discontinuity of T between zone 2 and zone 3 (interface L2)
% T2(L2)=(-q/2k2)*(L2^2)+B2*L2+E2
% T3(L2)=B3*L2+E3
% -k2*dT2/dx=hgap(T2(L2)-T3(L2))
% <==> -k2(-q/k2*L2+B2)=hgap(T2(L2)-T3(L2))
% <==> (k2+L2*hgap)B2+hgap*E2-L2*hgap*B3-hgap*E3=hgap*q/2k2*(L2^2)+q*L2
mat(4,1:6) =[0,0,k{2}(L(2))+L(2)*hgap,hgap,-L(2)*hgap,-hgap];
b(4) =hgap*(src{2}(L(2))/(2*k{2}(L(2))))*L(2)*L(2)+src{2}(L(2))*L(2);
% -k3*dT3/dx=hgap(T2(L2)-T3(L2))
% <==> -k3*B3=hgap(T2(L2)-T3(L2))
% <==> L2*hgap*B2+hgap*E2+(k3-L2*hgap)*B3-hgap*E3=hgap*q/2k2*(L2^2)
mat(5,1:6) =[0,0,L(2)*hgap,hgap,k{3}(L(2))-L(2)*hgap,-hgap];
b(5) =hgap*(src{2}(L(2))/(2*k{2}(L(2))))*L(2)*L(2);

a=mat\b';

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%