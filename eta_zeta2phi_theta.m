function [phi,theta] = eta_zeta2phi_theta(eta,zeta,r0,r1,r2)

r = cos(zeta)*(cos(eta)*r0+sin(eta)*r1)+sin(zeta)*r2;

phi = asin(r(3,1));
if r(2,1) < 0
    theta = 2*pi-acos(r(1,1)/cos(phi));
else
    theta = acos(r(1,1)/cos(phi));
end

