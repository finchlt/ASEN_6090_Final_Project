function [az,el,range] = ecef2azelrange2(r_sat,r_site,latgd,lon)

%==========================================================================
%==========================================================================
% [az,el,range] = ecef2azelrange2(r_sat,r_site,latgd,lon)
%
% Calculates the azimuth, elevation, and range of a satellite with respect
%  to an observation site.
%
%
% Author: Ben K. Bradley
% date: 11/15/2010
%
%
% INPUT:         Description                                         Units
%
%  r_sat      - position of satellite in ECEF frame                 [x y z]
%  r_site     - position of observing site in ECEF frame            [x y z]
%  latgd      - geodetic latitude of observation site          [-90,90] deg
%  lon        - longitude of observation site     [-180,180] or [0,360] deg
%
%
% OUTPUT:       
%    
%  az         - azimuth (degrees clockwise from North)          [0,360] deg
%  el         - elevation (degrees up from horizon)            [-90,90] deg
%  range      - distance from observation site to satellite                                    
%
%
% Coupling:
%
%  none
%
% References:
% 
%  [1] Vallado, D.A. "Fundamentals of Astrodynamics and Applications".
%         Third Edition. Microcosm Press. 2007.
%
%==========================================================================
%==========================================================================


%This algorithm comes from [1] p269

% Define Conversions ======================================================
halfpi  = pi*0.5;
rad2deg = 180/pi;
small   = 1e-8;

r_sat  = [r_sat(1); r_sat(2); r_sat(3)]; % make column vectors

r_site = [r_site(1); r_site(2); r_site(3)]; 


% Compute vector from observation site to satellite (ECEF) ================
rho_ecef = r_sat - r_site;  % [x; y; z]
range    = norm(rho_ecef);  % distance from observation site to satellite



% Convert ECEF to SEZ =====================================================
% % ROT2_lat = ROT2(halfpi - latgd/rad2deg);
% % ROT3_lon = ROT3(lon/rad2deg);
% % 
% % rho_sez  = ROT2_lat * ROT3_lon * rho_ecef;  % [S; E; Z]


cA = cosd(90-latgd);   sA = sind(90-latgd);
cB = cosd(lon);        sB = sind(lon);

T  = [cA*cB     cA*sB     -sA;
      -sB        cB        0;
      sA*cB     sA*sB      cA];
  
rho_sez = T * rho_ecef;  % [S; E; Z]  



% Compute Azimuth and Elevation ===========================================

% ELEVATION, rad ----------------------------------------------------

SE_mag = sqrt(rho_sez(1)*rho_sez(1) + rho_sez(2)*rho_sez(2));

if (SE_mag < small)  % directly over head
    
    el = sign(rho_sez(3)) * halfpi;   % +/- 90 degrees
    
else
    el = asin(rho_sez(3) / norm(rho_sez));
end

% AZIMUTH, rad ------------------------------------------------------

az = atan2(rho_sez(2)/SE_mag, -rho_sez(1)/SE_mag);



% Convert to degrees ======================================================
el = el * rad2deg;
az = az * rad2deg;

if (az < 0)
    az = az + 360;
end




