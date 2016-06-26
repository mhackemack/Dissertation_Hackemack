function F=hc_ss_x_DG_gap
% Solves the heat conduction equation in 1-D x-geometry using DFEM
% with T gap.
% An arbitrary number of material zones can be used but the analytical
% solution assumes 3 zones are used. The conductivities and the volumetric
% sources can be spatially dependent.

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
nel_zone = [ 4 10 2];

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
npar.ndofs = (npar.porder+1)*npar.nel;
% connectivity

gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=gn(iel-1,:) + npar.porder+1 ;
end
npar.gn=gn; clear gn;

% solve system
F = solve_fem3(dat,npar);
% plot
figure(1)

% create x values for plotting FEM solution
npar.xf=zeros(npar.nel*(npar.porder+1),1);
for iel=1:npar.nel
    ibeg = (iel-1)*(npar.porder+1)+1;
    iend = (iel  )*(npar.porder+1)  ;
    npar.xf(ibeg:iend)=linspace(npar.x(iel),npar.x(iel+1),npar.porder+1) ;
end

if length(nel_zone)==3
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

    plot(npar.xf,F,'.-',[x1 x2 x3],[y1 y2 y3],'r-'); hold all;
    title('1D steady-state heat conduction, 3 zones, with T gap, Cartesian coordinates')
    xlabel('Width (m)')
    ylabel('Temperature (C)')
    legend('FEM','Analytical','Location','northoutside','Orientation','horizontal')
else
    plot(npar.xf,F,'.-'); hold all;
    title('1D steady-state heat conduction, n zones, with T gap, Cartesian coordinates')
    xlabel('Width (m)')
    ylabel('Temperature (C)')
    legend('FEM','Location','northoutside','Orientation','horizontal')
end

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
g = npar.nelgap;
g1=g*(porder+1);
g2=g*(porder+1)+1;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+3)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*(porder+1);
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
% store shape set at edges
[be,dbedx] =feshpln([-1;1],porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      = int_domain (b rhs)

% loop over elements
for iel=1:nel
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
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end

xg=[1:(g-1) (g+1):(nel-1)] ;

% loop on interior edges
for ie=xg
    % left (interior)/right (exterior) elements
	ieli = ie;
	iele = ie+1;
    % values for 2 cells 1= left of edge, 2= right of edge
    ce1=npar.iel2zon(ieli);
    ce2=npar.iel2zon(iele);
    % element extremities
    x0=npar.x(ie);
    x1=npar.x(ie+1);
    x2=npar.x(ie+2);
    % element lengths
    d1=x1-x0;
    d2=x2-x1;
    % conductivity values
    k1=dat.k{ce1}(x1);
    k2=dat.k{ce2}(x1);
    % compute local matrices + load vector
    mee=(k2/d2)*dbedx(1,:)'*be(1,:);
    mei=(k1/d1)*dbedx(end,:)'*be(1,:);
    mie=(-k2/d2)*dbedx(1,:)'*be(end,:);
    mii=(-k1/d1)*dbedx(end,:)'*be(end,:);
    kap=2*(k1/d1+k2/d2);
    % assemble
    A(gn(ieli,:),gn(ieli,:)) = A(gn(ieli,:),gn(ieli,:)) + mii' + mii;
    A(gn(ieli,:),gn(iele,:)) = A(gn(ieli,:),gn(iele,:)) + mie' + mei;
    A(gn(iele,:),gn(ieli,:)) = A(gn(iele,:),gn(ieli,:)) + mei' + mie;
    A(gn(iele,:),gn(iele,:)) = A(gn(iele,:),gn(iele,:)) + mee' + mee;
    % assemble
    A(gn(ieli,end),gn(ieli,end)) = A(gn(ieli,end),gn(ieli,end))+kap;
    A(gn(ieli,end),gn(iele,1))   = A(gn(ieli,end),gn(iele,1))  -kap;
    A(gn(iele,1),gn(ieli,end))   = A(gn(iele,1),gn(ieli,end))  -kap;
    A(gn(iele,1),gn(iele,1))     = A(gn(iele,1),gn(iele,1))    +kap;
end

A(g1,g1)=A(g1,g1)+dat.hgap;
A(g1,g2)=A(g1,g2)-dat.hgap;
A(g2,g1)=A(g2,g1)-dat.hgap;
A(g2,g2)=A(g2,g2)+dat.hgap;

% apply natural BC
% element extremities
x0=npar.x(1);
x1=npar.x(2);
xnm=npar.x(end-1);
xn=npar.x(end);
% element lengths
d1=x1-x0;
dn=xn-xnm;
% conductivity values
k1=dat.k{npar.iel2zon(1)}(x1);
kn=dat.k{npar.iel2zon(end)}(xn);
% compute load vector
kap1=8*k1/d1;
kapn=8*kn/dn;
    
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+dat.hcv;
        rhs(1)=rhs(1)+dat.hcv*dat.bc.left.C;
case 2 % Dirichlet
        minus_one=1;
        A(gn(1,:),gn(1,:))=A(gn(1,:),gn(1,:))+kap1*be(minus_one,:)'*be(minus_one,:)...
            +2*k1*be(minus_one,:)'   *dbedx(minus_one,:)/d1...
            +2*k1*dbedx(minus_one,:)'*be(minus_one,:)   /d1;
        rhs(gn(1,:))=rhs(gn(1,:))+dat.bc.left.C*...
            ( kap1*be(minus_one,:)' +2*k1*dbedx(minus_one,:)'/d1 );
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+dat.hcv;
        rhs(n)=rhs(n)+dat.hcv*dat.bc.rite.C;
    case 2 % Dirichlet
        plus_one =2;
        A(gn(end,:),gn(end,:))=A(gn(end,:),gn(end,:))...
            +kapn*be(plus_one,:)'*be(plus_one,:)...
            -2*kn*be(plus_one,:)'   *dbedx(plus_one,:)/dn...
            -2*kn*dbedx(plus_one,:)'*be(plus_one,:)   /dn;
        rhs(gn(end,:))=rhs(gn(end,:))+dat.bc.rite.C*...
            ( kapn*be(plus_one,:)' -2*kn*dbedx(plus_one,:)'/dn );
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