function [mp1] = MP1(P1,L1,L2)

% some constants
c = 299792458; %speed of light in m/s
omegaE = 7.2921151467e-5; % 2*pi/86400; Earth rotation Rate
f_L1 = 1.57542e3;         % L1 freq in MHz
f_L2 = 1.22760e3;         % L2 freq in MHz
lambda_L1 = c/(f_L1*1e6); % L1 wavelength in m
lambda_L2 = c/(f_L2*1e6); % L2 wavelength in m

mp1 = P1 - ( ((f_L1^2+f_L2^2)/(f_L1^2-f_L2^2))*L1*lambda_L1) + ( (2*f_L2^2/(f_L1^2-f_L2^2))*L2*lambda_L2 );

end