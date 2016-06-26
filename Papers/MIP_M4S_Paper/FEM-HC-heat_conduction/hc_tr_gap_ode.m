clear all; close all; clc;

dat.k{1}=@k_Zr;
dat.k{2}=@k_fuel;
dat.k{3}=@k_clad;

dat.hgap= 15764    /1e4; % units: W/cm^2-C
dat.hcv = 2574     /1e4; % units: W/cm^2-C
dat.width=[0.003175 0.0174115 0.0179195]*100; % converting in cm

bc.rite.type=1;   % 0=neumann, 1=robin, 2=dirichlet
bc.rite.Tcool=49.51; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)

% from Chance:
node_height= 38.1/15; % units: cm
q_prime = 1099.2/node_height; % units: W/cm
q_prime = 438.72; % based on Tclad outer - Tcool from Fig.8 in Chance's)

%  check
if bc.rite.type ~=1
    error('analytical solution only coded for Robin');
end

% shortcuts. 1 is a dummy radial position since properties are uniform
kz = dat.k{1}(1)/100; % units: W/cm-C
kf = dat.k{2}(1)/100; % units: W/cm-C
kc = dat.k{3}(1)/100; % units: W/cm-C
Tcool = bc.rite.Tcool;
R1=dat.width(1);
R2=dat.width(2);
R3=dat.width(3);
hgap =dat.hgap;
hconv=dat.hcv;

%%%%%%%%%%%%%%%%
%%%% clad
%%%%%%%%%%%%%%%%
T_outer_clad = q_prime/(2*pi*R3*hconv) + Tcool

rc=linspace(R2,R3,3);
Tc = q_prime/(2*pi*kc) * log(R3./rc) + T_outer_clad;

%%%%%%%%%%%%%%%%
%%%% fuel
%%%%%%%%%%%%%%%%
Tc(1)
T_outer_fuel = q_prime/(2*pi*R2*hgap) + Tc(1)

rf=linspace(R1,R2,10);
Tf = q_prime / ( 2*pi*kf*(R2^2-R1^2) ) * ( (R2^2-rf.^2)/2  - R1^2*log(R2./rf) )  + T_outer_fuel;

%%%%%%%%%%%%%%%%
%%%% zr rod
%%%%%%%%%%%%%%%%
Tf(1)
rz=linspace(0,R1,3);
Tz = ones(1,length(rz))*Tf(1);

%%%%%%%%%%%%%%%%
%%%% single vector

r=[rz rf rc R3*1.01 R3*1.1];
T=[Tz Tf Tc Tcool Tcool];


plot(r,T,'.-');

%%%%%%%%%%%%%%%%%
% thermal resistance: Rth Delta_t = q'
Rth = q_prime / ( Tf(1) - Tcool );

% 0-D time-dependent
% rhocp dT/dt + Rth(T-Tcool)=q
% rhocp (Tnew-Told)/dt + Rth(Tnew-Tcool)=q
rhocp=cp_fuel(1)*rho_fuel(1)/1e6; % kg/m3 . J/kg/C * conversion(m/cm)^3 = J/cm3-C 
rhocpf=rhocp*pi*(R2^2-R1^2); % time surface of fuel, gives units of J/cm-C
rhocp =rhocpf + cp_clad(1)*rho_clad(1)/1e6*pi*(R3^2-R2^2 + R1^2-0^2); % rhocp=rhocp*pi*(R3^2); % time surface of fuel, gives units of J/cm-C
% rhocp=rhocpf;
% rhocp = cp_fuel(1)*rho_fuel(1)/1e6*pi*R3^2;
% [rhocp*pi*(R2^2-R1^2)  cp_clad(1)*rho_clad(1)/1e6*pi*(R3^2-R2^2 + R1^2-0^2)]


time_end = 1e2;
ntimes=50;
dt = time_end/ntimes;
Told=Tcool; % initial value
T=zeros(ntimes+1,1);T(1)=Told;
for it=1:ntimes
    Tnew = ( rhocp/dt*Told + q_prime + Rth*Tcool ) / ( rhocp/dt+Rth );
    T(it+1)=Tnew;
    Told=Tnew;
end

t=linspace(0,time_end,ntimes+1);
figure(20);
plot(t,T);
xlabel('time, s');
ylabel('T, C');
