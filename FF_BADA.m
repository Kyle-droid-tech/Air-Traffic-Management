function [delt, FF] = FF_BADA(x1, x2, h1, h2, TAS1, TAS2,  OPF)

ft2m  = 0.3048;
kt2ms = 0.5144;
g0 = 9.80665;

Swing = OPF(1,1);
mass = 190*10^3;  % OPF(2,1);
CD0    = OPF(10,1);
CD2    = OPF(11,1);
Cf1    = OPF(12,1);
Cf2    = OPF(13,1);
Cf3    = OPF(14,1);
Cf4    = OPF(15,1);
% Cfcr   = OPF(16,1);
% CTc1   = OPF(17,1);
% CTc2   = OPF(18,1);
% CTc3   = OPF(19,1);
% CTc4   = OPF(20,1);
% CTc5   = OPF(21,1);
k2=1.0/60.0/10.0^3; 
k3=1.0/60.0;
f0 = k3*Cf3;
f1 = k3*Cf3/Cf4;
f2 = k2*Cf1;
f3 = k2*Cf1/Cf2;

dh = h2-h1;
h_m = (h1+h2)/2;
[~,rho_m,~,~] = ISA(h_m*ft2m);
dx = x2-x1;

% [~,rho1,p1,~] = ISA(h1*ft2m);
% [~,rho2,p2,~] = ISA(h2*ft2m);
% TAS1 = cas2tas(rho1,p1*100,CAS1*kt2ms);
% TAS2 = cas2tas(rho2,p2*100,CAS2*kt2ms);

dTAS = TAS2-TAS1;
TAS_m = (TAS1+TAS2)/2; 
gam = atan2(dh*ft2m, dx);

delt = dx/TAS_m; % No wind condition
CL = mass*g0*cos(gam)/(1/2*rho_m*TAS_m^2*Swing);
CD = CD0+CD2*CL^2;
D  = 1/2*rho_m*TAS_m^2*Swing*CD;
T  = mass*dTAS/delt+D+mass*g0*sin(gam);
FFnom = (f2+f3*TAS_m/kt2ms)*T;
FFmin = f0-f1*h_m/ft2m;

FF = max(FFnom, FFmin);
