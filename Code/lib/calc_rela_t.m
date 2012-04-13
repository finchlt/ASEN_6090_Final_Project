function [rela_t] = calc_rela_t(curreph,t)
mu = 3986004.418e8; % m^3/s^2 earth's gravitational const
c = 299792458; %speed of light in m/s

a = (curreph(5))^2;
e = curreph(4);
meanAnom = curreph(2);
toe = curreph(17);

% Compute SV mean motion n and mean anomaly M.
n=sqrt(mu./(a.^3));
M=meanAnom+n.*(t-toe); 
M=mod(M, 2*pi);

% Compute SV orbit eccentricity anomaly E iteratively. 
err=1; 
E=M; 
while (err>1e-5)
oldE=E; 
E=M+(e).*sin(E); 
err=max(abs(E-oldE));
end

rela_t = (-2/c^2)*sqrt(mu*a)*(e*sin(E));
end