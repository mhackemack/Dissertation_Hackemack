function F=hc_tr_r_DG
% Solves the time-dependent heat conduction equation in 1-D r-geometry 
% using DFEM without T gap.
% The conductivities and the volumetric sources can be spatially dependent.

% clear the console screen
clc; clear all; close all;
% load the data structure with info pertaining to the physical problem
dat.k{1}=@k_Zr; % W/m-K
dat.k{2}=@k_fuel; % W/m-K
dat.k{3}=@k_clad; % W/m-K
dat.esrc{1}=@zero_functiont; % W/m3
dat.esrc{2}=@esrct; % W/m3
dat.esrc{3}=@zero_functiont; % W/m3
dat.rho{1}=@rho_Zr; % kg/m3
dat.rho{2}=@rho_fuel; % kg/m3
dat.rho{3}=@rho_clad; % kg/m3
dat.cp{1}=@cp_Zr; % J/kg-K
dat.cp{2}=@cp_fuel; % J/kg-K
dat.cp{3}=@cp_clad; % J/kg-K

dat.hcv=1612.414; % W/m^2-C
dat.hgap=15764; % W/m^2-C
dat.width=[0.003175 0.0174115 0.0179195]; % m
dat.duration=10000; % s
dat.Tinit=30; % Celsius
bc.rite.type=1; % 0=neumann, 1=robin, 2=dirichlet
bc.rite.C=50; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
dat.bc=bc; clear bc;

% number of time points
npar.n_time_steps=100; % total number of time steps to run
npar.delta_t=dat.duration/npar.n_time_steps;

% number of the curves to plot, max 14 curves
npar.curve=10; 

nel_zone = [ 20 100 5];

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
T = time_solve(dat,npar);
% plot
plot_solution(dat,npar,T)

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T=time_solve(dat,npar)

% shortcuts
nel   = npar.nel;
porder= npar.porder;
n_time_steps=npar.n_time_steps;
delta_t=npar.delta_t;
curve=npar.curve;

n=nel*(porder+1);

% initial condition
T=zeros(n,curve);
T_old=dat.Tinit*ones(n,1);

nummax=n_time_steps+1;
increment=round(nummax/curve);

% loop over time
end_time=0;
for it=1:curve
    end_time=end_time+increment*delta_t;
    % solve the FME system for a given time step
    T_old=assemble_solve(dat,npar,end_time,T_old);
	% save solutions at different time points
    T(:,it)=T_old;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T=assemble_solve(dat,npar,time,T_old)
% assemble the matrix, the rhs, apply BC and solve

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
delta_t=npar.delta_t;
L = dat.width;

% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+3)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*(porder+1);
% allocate memory
A=spalloc(n,n,nnz);
rhs=zeros(n,1);
T=zeros(n,1);

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
    Mival=(x1+x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    my_zone=npar.iel2zon(iel);
    d=dat.k{my_zone}(x);
    q=dat.esrc{my_zone}(x,time);
    rho=dat.rho{my_zone}(x);
    cp=dat.cp{my_zone}(x);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            m(i,j)= dot((Mival.*rho.*cp.*wq.*b(:,i)+Jac.*xq.*rho.*cp.*wq.*b(:,i)), b(:,j));
            k(i,j)= dot((Mival.*d.*wq.*dbdx(:,i)+Jac.*xq.*d.*wq.*dbdx(:,i)), dbdx(:,j));
        end
        f(i)= dot((Mival.*q.*wq+Jac.*xq.*q.*wq), b(:,i));
    end
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + m*Jac/delta_t + k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + m*Jac*T_old(gn(iel,:))/delta_t + f*Jac;
end

% loop on interior edges
for ie=1:(npar.nel-1)
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
    mee=x1*(k2/d2)*dbedx(1,:)'*be(1,:);
    mei=x1*(k1/d1)*dbedx(end,:)'*be(1,:);
    mie=x1*(-k2/d2)*dbedx(1,:)'*be(end,:);
    mii=x1*(-k1/d1)*dbedx(end,:)'*be(end,:);
    kap=2*x1*(k1/d1+k2/d2);
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

% apply natural BC
% element extremities
xnm=npar.x(end-1);
xn=npar.x(end);
% element lengths
dn=xn-xnm;
% conductivity values
kn=dat.k{npar.iel2zon(end)}(xn);
% compute load vector
kapn=8*kn/dn;

switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+L(3)*dat.hcv;
        rhs(n)=rhs(n)+L(3)*dat.hcv*dat.bc.rite.C;
    case 2 % Dirichlet
        plus_one =2;
        A(gn(end,:),gn(end,:))=A(gn(end,:),gn(end,:))...
            +L(2)*(kapn*be(plus_one,:)'*be(plus_one,:)...
            -2*kn*be(plus_one,:)'   *dbedx(plus_one,:)/dn...
            -2*kn*dbedx(plus_one,:)'*be(plus_one,:)   /dn);
        rhs(gn(end,:))=rhs(gn(end,:))+L(2)*dat.bc.rite.C*...
            ( kapn*be(plus_one,:)' -2*kn*dbedx(plus_one,:)'/dn );
end

% solve
T=A\rhs;

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
function plot_solution(dat,npar,T)
% plot solution for different time points            

figure(1); hold all;

% create x values for plotting FEM solution
npar.xf=zeros(npar.nel*(npar.porder+1),1);
for iel=1:npar.nel
    ibeg = (iel-1)*(npar.porder+1)+1;
    iend = (iel  )*(npar.porder+1)  ;
    npar.xf(ibeg:iend)=linspace(npar.x(iel),npar.x(iel+1),npar.porder+1) ;
end

nummax=npar.n_time_steps+1;
increment=round(nummax/npar.curve);
color=1;
end_time=0;
for num=1:npar.curve
	end_time=end_time+increment*npar.delta_t;
	c=val_color(color);
	plot(npar.xf,T(:,num),c);
	name(color)=cellstr(strcat('time=',num2str(end_time)));
	legend_graph(color,:)=name(color);
	color=color+1;
	legend(legend_graph,'Location','southwest');
end

title('1D transient heat conduction, without T gap, cylindrical coordinates')
xlabel('Width (m)')
ylabel('Temperature (C)')
grid on

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%