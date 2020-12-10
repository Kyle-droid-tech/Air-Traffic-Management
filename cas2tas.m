function vtas = cas2tas(rho,p,vcas)

%% CAS to TAS transfer function for 4D optimization

p0=101325; %Pa
rho0=1.225; %kg/m^3
kappa = 1.4;
mu = (kappa-1)/kappa;

c1   = (1+mu*rho0*vcas^2/(2*p0))^(1/mu)-1;
c2   = (1+(p0/p)*c1)^mu-1;
vtas = sqrt(2*p*c2/(mu*rho));
