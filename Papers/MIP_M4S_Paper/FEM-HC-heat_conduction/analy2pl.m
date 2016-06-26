function F=analy2pl
% Solves analytically the heat conduction equation in 1-D r-geometry
% with T gap.

% clear the console screen
clc; close all;
% load the data structure with info pertaining to the physical problem
dat.k{1}=@(r) 18;
dat.k{2}=@(r) 16;
dat.k{3}=@(r) 16;
dat.esrc{1}=@(r) 0;
dat.esrc{2}=@(r) 5000000;
dat.esrc{3}=@(r) 0;
dat.hgap=15764;
dat.hcv=20000;
dat.width=[0.006 0.034823 0.039];
bc.rite.type=1;
bc.rite.C=400;
dat.bc=bc; clear bc;

verif_hc_eq(dat);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function verif_hc_eq(dat)

k=dat.k; src=dat.esrc; hgap=dat.hgap; hcv=dat.hcv; L=dat.width;

% general form of the solution:
% Zone 1 : T1 = E1
% dT1/dr= 0
% Zone 2 : T2 = -q/(4*k2)*(r^2) + B2*ln(r) + E2
% dT2/dr= -q/2k2*r + B2/r
% Zone 3 : T3 = B3*ln(r) + E3
% dT3/dr= B3/r

switch dat.bc.rite.type
    case 0 % Neumann
        % k3*du3/dn=C on the right becomes: k3*du3/dx=C
        % <==> k3*B3/L3=C <==> B3=C*L3/k3
        mat(5,1:5) =[0,0,0,1,0];
        b(5) = dat.bc.rite.C*L(3) / k{3}(L(3));
    case 1 % Robin
        % u3+k3/hcv*du3/dn =C on the right becomes: u3+k3/hcv*du3/dx =C
        % <==> (ln(L3)+k3/(hcv*L3))B3+E3=C
        mat(5,1:5) =[0,0,0,log(L(3))+k{3}(L(3))/(hcv*L(3)),1];
        b(5) = dat.bc.rite.C;
    case 2 % Dirichlet
        % u3=C <==> B3*ln(L3)+E3=C
        mat(5,1:5) =[0,0,0,log(L(3)),1];
        b(5) = dat.bc.rite.C;
end

% fixed conditions
% continuity of T and flux between zone 1 and zone 2 (interface L1)
% T1(L1)=T2(L1) <==> B1=(-q/4k2)*(L1^2)+B2*ln(L1)+E2
% <==> -E1+B2*ln(L1)+E2=(q/4k2)*(L1^2)
mat(1,1:5) =[-1,log(L(1)),1,0,0];
b(1) =src{2}(L(1))/(4*k{2}(L(1)))*L(1)*L(1);
% phi1(L1)=phi2(L1) <==> 0=k2*((-q/2k2)*L1+B2/L1)
% <==> B2=(q/2k2)*L1*L1
mat(2,1:5) =[0,1,0,0,0];
b(2) =src{2}(L(1))/(2*k{2}(L(1)))*L(1)*L(1);

% discontinuity of T between zone 2 and zone 3 (interface L2)
% T2(L2)=(-q/4k2)*(L2^2)+B2*ln(L2)+E2
% T3(L2)=B3*ln(L2)+E3
% -k2*dT2/dr=hgap(T2(L2)-T3(L2))
% <==> -k2(-q/2k2*L2+B2/L2)=hgap(T2(L2)-T3(L2))
% <==> (k2/(hgap*L2)+ln(L2))*B2+E2-ln(L2)*B3-E3=q/4k2*(L2^2)+q*L2/(2*hgap)
mat(3,1:5) =[0,k{2}(L(2))/(hgap*L(2))+log(L(2)),1,-log(L(2)),-1];
b(3) =(src{2}(L(2))/(4*k{2}(L(2))))*L(2)*L(2)+src{2}(L(2))*L(2)/(2*hgap);
% -k3*dT3/dr=hgap(T2(L2)-T3(L2))
% <==> -k3*B3/L2=hgap(T2(L2)-T3(L2))
% <==> ln(L2)*B2+E2+(k3/(L2*hgap)-ln(L2))*B3-E3=q/4k2*(L2^2)
mat(4,1:5) =[0,log(L(2)),1,k{3}(L(2))/(L(2)*hgap)-log(L(2)),-1];
b(4) =(src{2}(L(2))/(4*k{2}(L(2))))*L(2)*L(2);

% get coefficient for the analytical solution
a=mat\b';
r1=linspace(0,L(1));
r2=linspace(L(1),L(2));
r3=linspace(L(2),L(3));
y1=a(1)+0*r1;
y2=-src{2}(r2)/(4*k{2}(r2))*(r2.^2)+a(2)*log(r2)+a(3);
y3=a(4)*log(r3)+a(5);

plot(r1,y1,r2,y2,r3,y3); hold all;
title('1D heat conduction problem')
xlabel('Width')
ylabel('Temperature')

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%