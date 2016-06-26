function F=analy2
% Solves analytically the heat conduction equation in 1-D x-geometry
% with T gap.

% clear the console screen
clc; close all;
% load the data structure with info pertaining to the physical problem
global dat
dat.k{1}=@(x) 18;
dat.k{2}=@(x) 16;
dat.k{3}=@(x) 16;
dat.esrc=5000000;
dat.hgap=15764;
dat.hcv=20000;
dat.width=[0.003175 0.034823 0.036];
bc.left.type=1; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=400; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
bc.rite.type=1;
bc.rite.C=80;
dat.bc=bc; clear bc;

verif_hc_eq;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function verif_hc_eq
global dat

k1=dat.k{1}; 
k2=dat.k{2}; 
k3=dat.k{3}; 
src=dat.esrc;
hgap=dat.hgap; hcv=dat.hcv; L=dat.width; 

% general form of the solution:
% Zone 1 : T1 = B1*x + E1
% dT1/dx= B1
% Zone 2 : T2 = -q/(2*k2)*(x.^2) + B2*x + E2 = Y2*(x.^2) + B2*x + E2
% dT2/dx= -q/k2*x + B2 = Y1*x + B2
% Zone 3 : T3 = B3*x + E3
% dT3/dx= B3

% definition of 2 coefficients
Y1=-src/k2(0);
Y2=-src/(2*k2(0));

switch dat.bc.left.type
    case 0 % Neumann
        % k1*du1/dn=C on the left becomes: -k1*du1/dx=C 
        % <==> -k1*B1=C <==> B1=-C/k1
        mat(1,1:6) =[1,0,0,0,0,0];
        b(1) = -dat.bc.left.C / k1(0);
    case 1 % Robin
        % u1+k1/hcv*du1/dn =C on the left becomes: u1-k1/hcv*du1/dx =C
        % <==> -k1/hcv*B1+E1=C
        mat(1,1:6) =[-k1(0)/hcv,1,0,0,0,0];
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
        b(6) = dat.bc.rite.C / k3(0);
    case 1 % Robin
        % u3+k3/hcv*du3/dn =C on the right becomes: u3+k3/hcv*du3/dx =C
        % <==> (L3+k3/hcv)B3+E3=C
        mat(6,1:6) =[0,0,0,0,L(3)+k3(0)/hcv,1];
        b(6) = dat.bc.rite.C;
    case 2 % Dirichlet
        % u3=C <==> B3*L+E3=C
        mat(6,1:6) =[0,0,0,0,L(3),1];
        b(6) = dat.bc.rite.C;
end

% fixed conditions
% continuity of T and flux between zone 1 and zone 2 (interface L1)
% T1(L1)=T2(L1) <==> B1*L1+E1=Y2*(L1^2)+B2*L1+E2
% <==> B1*L1+E1-B2*L1-E2=Y2*(L1^2)
mat(2,1:6) =[L(1),1,-L(1),-1,0,0];
b(2) =Y2*L(1)*L(1);
% phi1(L1)=phi2(L1) <==> k1*B1=k2*(Y1*L1+B2)
% <==> (k1/k2)*B1-B2=Y1*L1
mat(3,1:6) =[k1(0)/k2(0),0,-1,0,0,0];
b(3) =Y1*L(1);

% discontinuity of T between zone 2 and zone 3 (interface L2)
% T2(L2)=Y2*(L2^2)+B2*L2+E2
% T3(L2)=B3*L2+E3
% -k2*dT2/dx=hgap(T2(L2)-T3(L2))
% <==> -k2(-q/k2*L2+B2)=hgap(T2(L2)-T3(L2))
% <==> (k2+L2*hgap)B2+hgap*E2-L2*hgap*B3-hgap*E3=-hgap*Y2*(L2^2)+q*L2
mat(4,1:6) =[0,0,k2(0)+L(2)*hgap,hgap,-L(2)*hgap,-hgap];
b(4) =-hgap*Y2*L(2)*L(2)+src*L(2);
% -k3*dT3/dx=hgap(T2(L2)-T3(L2))
% <==> -k3*B3=hgap(T2(L2)-T3(L2))
% <==> L2*hgap*B2+hgap*E2+(k3-L2*hgap)*B3-hgap*E3=-hgap*Y2*(L2^2)
mat(5,1:6) =[0,0,L(2)*hgap,hgap,k3(0)-L(2)*hgap,-hgap];
b(5) =-hgap*Y2*L(2)*L(2);

% get coefficient for the analytical solution
a=mat\b';
x1=linspace(0,L(1));
x2=linspace(L(1),L(2));
x3=linspace(L(2),L(3));
y1=a(1)*x1+a(2);
y2=Y2*(x2.^2)+a(3)*x2+a(4);
y3=a(5)*x3+a(6);

plot(x1,y1,x2,y2,x3,y3); hold all;
title('1D heat conduction problem')
xlabel('Width')
ylabel('Temperature')

return
end