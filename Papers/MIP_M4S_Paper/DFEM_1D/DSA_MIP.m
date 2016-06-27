function A=DSA_MIP()
% global strucures
global npar mats dat

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=nel*(porder+1);
% Build volumetric terms
% ------------------------------------------------------------------------------
A = mats.elem.k + mats.elem.ma;
% Build interior edge terms
% ------------------------------------------------------------------------------
for ie=1:(npar.nel-1)
    % left (interior)/right (exterior) elements
    ieli = ie;
    iele = ie+1;
    % values for 2 cells 1= left of edge, 2= right of edge
    ze1=npar.iel2zon(ieli);
    ze2=npar.iel2zon(iele);
    % element extremities
    x0=npar.x(ie);
    x1=npar.x(ie+1);
    x2=npar.x(ie+2);
    % element lengths
    dx1=x1-x0;
    dx2=x2-x1;
    % conductivity values
    D1=dat.D(ze1);
    D2=dat.D(ze2);
    % compute local matrices
    gL  = mats.elem.egr{1}{iele};
    gR  = mats.elem.egr{2}{ieli};
    gCL = mats.elem.egrC{1}{iele};
    gCR = mats.elem.egrC{2}{ieli};
    
    mee =  (D2/2)*gL;
    mei =  (D1/2)*gCR;
    mie = (-D2/2)*gCL;
    mii = (-D1/2)*gR;
%     mee = (D2/2)*dbedx(1,:)'*be(1,:);
%     mei = (D1/2)*dbedx(end,:)'*be(1,:);
%     mie = (-D2/2)*dbedx(1,:)'*be(end,:);
%     mii = (-D1/2)*dbedx(end,:)'*be(end,:);
    kap=get_kappa(4,porder,[D1,D2],[dx1,dx2],0);
%     kap=2*(D1/d1+D2/d2);
    % assemble gradient terms
    A(gn(ieli,:),gn(ieli,:)) = A(gn(ieli,:),gn(ieli,:)) + (mii' + mii);
    A(gn(ieli,:),gn(iele,:)) = A(gn(ieli,:),gn(iele,:)) + (mie' + mei);
    A(gn(iele,:),gn(ieli,:)) = A(gn(iele,:),gn(ieli,:)) + (mei' + mie);
    A(gn(iele,:),gn(iele,:)) = A(gn(iele,:),gn(iele,:)) + (mee' + mee);
    % assemble mass terms
    A(gn(ieli,end),gn(ieli,end)) = A(gn(ieli,end),gn(ieli,end))+kap;
    A(gn(ieli,end),gn(iele,1))   = A(gn(ieli,end),gn(iele,1))  -kap;
    A(gn(iele,1),gn(ieli,end))   = A(gn(iele,1),gn(ieli,end))  -kap;
    A(gn(iele,1),gn(iele,1))     = A(gn(iele,1),gn(iele,1))    +kap;
end

% Left Boundary Condition
% ------------------------------------------------------------------------------
gn1 = gn(1,:); gn2 = gn(2,:);
zn1 = npar.iel2zon(1);
zn2 = npar.iel2zon(2);
D1  = dat.D(zn1);
D2  = dat.D(zn2);
% x0  = npar.x(1);
% x1  = npar.x(2);
% x2  = npar.x(3);
% dx1 = x1-x0;
% dx2 = x2-x1;
% kap = get_kappa(4,porder,D,dx1,1);
% mb  = (D1/2)*mats.elem.egr{1}{1};
nL  = -1; nR  = 1;
gL  = (D1/2)*mats.elem.egr{1}{1};
gR  = (D1/2)*mats.elem.egr{2}{1};
% original dirichlet
% A(1,1) = A(1,1) + kap;
% A(gn1,gn1) = A(gn1,gn1) - n*(mb + mb');
% modified condition
A(1,1) = A(1,1) + 1/2;
A(gn1,gn1) = A(gn1,gn1) + nL*gL;
A(gn1,gn1) = A(gn1,gn1) + nR*gR;

% Right Boundary Condition
% ------------------------------------------------------------------------------
gn1 = gn(end-1,:); gn2 = gn(end,:);
zn1 = npar.iel2zon(end-1);
zn2 = npar.iel2zon(end);
D1  = dat.D(zn1);
D2  = dat.D(zn2);
% x0  = npar.x(end-2);
% x1  = npar.x(end-1);
% x2  = npar.x(end);
% dx1 = x1-x0;
% dx2 = x2-x1;
% kap = get_kappa(4,porder,D,dx2,1);
% mb  = (D/2)*mats.elem.egr{2}{end};
nL  = 1; nR  = 1;
gL  = (D1/2)*mats.elem.egr{1}{end};
gR  = (D1/2)*mats.elem.egr{2}{end};
% original dirichlet
% A(end,end) = A(end,end) + kap;
% A(gnn,gnn) = A(gnn,gnn) - n*(mb + mb');
% modified condition
A(end,end) = A(end,end) + 1/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_kappa(C,p,D,h,eflag)
c = C*(1+p)*p;
if eflag == 0
    out = c/2*(D(1)/h(1) + D(2)/h(2));
else
    out = c*D/h;
end
out = max(out, 0.25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
