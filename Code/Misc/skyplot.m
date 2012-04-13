function [h,theta,rho] = skyplot(elev,azi,prn,h,r_or_d)
% SKYPLOT, plot of satellite path
%
%   [h] = SKYPLOT(elev,azi,prn,h,r_or_d)
%   Generates a skyplot of a satellite path given its azimuth and
%   elevation
%   r_or_d, specifies the angles in radians or degree (1->rad,0->deg;default degree)
%   prn is the SV name (string)
%   h is the figure window handle

if size(elev) ~= size(azi)
    disp('Invalid Elevation/Azimuth values')
    return
end

make_bases = 0;

if nargin == 2
    r_or_d = 0;
    
    prn = ' ';
    
    h = figure;
    make_bases = 1;
elseif nargin == 3
    r_or_d = 0;
    
    h = figure;
    make_bases = 1;
elseif nargin == 4
    r_or_d = 0;
end


if(make_bases == 1)
    theta = 0:0.001:2*pi;
    figure(h), hold on
    for r = 0.2:0.2:1
        figure(h),plot(r.*cos(theta),r.*sin(theta),'k--');
        %e_val = r*(-pi/2)+pi/2;
        %text(r*cos(theta(1)),r*sin(theta(1)),num2str(e_val))
    end
    figure(h), hold off
    figure(h), axis off
end

if r_or_d == 0
    elev = elev*pi/180;
    azi = azi*pi/180;
end

theta = -azi + pi/2;
rho = (elev-pi/2)*(-1)/(pi/2);

figure(h),hold on,
figure(h),plot(rho.*cos(theta),rho.*sin(theta))
text(rho(10)*cos(theta(10)),rho(10)*sin(theta(10)),prn);
plot(rho(10)*cos(theta(10)),rho(10)*sin(theta(10)),'go','MarkerSize',25);
figure(h),hold off