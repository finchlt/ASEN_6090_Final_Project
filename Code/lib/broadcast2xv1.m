function [health,x,v,curreph] = broadcast2xv1(ephem_all,t_input,prn)

%==========================================================================
%==========================================================================
% [health,x,v] = broadcast2xv(ephem_all,t_input,prn)
%
% Calculates the position and velocity of a GPS satellite 
%  from an ephemeris matrix (see read_GPSbroadcast.m).  
%
%
% Author: Ben K. Bradley
% date: 07/19/2009
% Additional features in the original version deleted P. Axelrad
% date: 10/17/2011
%
%
% INPUT:               Description                                  Units
%
%  ephem_all    - matrix of gps satellite orbit parameters
%
%        = [prn M0 delta_n e sqrt_a Loa i perigee ra_rate i_rate Cuc Cus...
%               Crc Crs Cic Cis Toe IODE GPS_week Toc Af0 Af1 Af2 0 health]
%
%
%  t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
%  prn          - PRN to compute values for (one satellite only)                       
%
%
% OUTPUT:       
%    
%  health       - health of satellite (0=good)                     (nx1)
%  x            - position of satellite (ECEF)         [x y z]   m (nx3)
%  v            - velocity of satellite (ECEF)         [vx vy vz] m/s (nx3)
%  accel        - acceleration of satellite (ECEF)     [ax ay az] m/s^2 (nx3)
%                                     
%
% Coupling:
%
%   mean2eccentric.m
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
%         Determination using the Broadcast Ephemeris". The Journal of
%         Navigation. (2006), 59, 293-305.
%            < http://journals.cambridge.org/action/displayAbstract;jsess ...
%                ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
%
%   [3] skyplot.cpp by the National Geodetic Survey
%          < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >
%
%==========================================================================
%==========================================================================


%NOTE: Numbered equations in the code (eg.  Eq. 21) correspond to equations
%        in the [2] reference.


% Load GPS Accepted WGS-84 Constants ======================================
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
wE  = 7.2921151467e-5; % WGS-84 value, rad/s 
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s
%PI = 3.1415926535898; % accepted GPS value for pi



% Initialize Output Variables for Speed ===================================
sz = size(t_input,1);

x          = ones(sz,3) * NaN;
v          = x; 
accel      = x;
health     = ones(sz,1) * NaN; 
satClkCorr = health;
relcorr    = health;

t_sec = t_input(:,2); 
gpswk = t_input(:,1);                                          



% Start Main Calculation Loop =============================================
%==========================================================================

aaa  = find(ephem_all(:,1) == prn);  % aaa is vector containing row numbers of ephem_all that are for sat.no. 'index' 
sat_ephem = ephem_all(aaa,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'


curreph = zeros(length(t_sec),size(ephem_all,2));
for ttt = 1:length(t_sec) %loop through all times

%     rrr = find(gpswk(ttt) == sat_ephem(:,19));
%     sat_ephem = sat_ephem(rrr,:);

%     jjj = max( find(t_sec(ttt) >= sat_ephem(:,17)) ); % sat_ephem(:,17) = toe (sec into GPS week) of each entry
                                                     % jjj = row of specific sat. ephem. data with epoch closest to input time
    [v jjj] = min( abs(t_sec(ttt) - sat_ephem(:,17)) );
    if isempty(jjj),break,end 
    
    curreph(ttt,:) = sat_ephem(jjj,:);

    toe = sat_ephem(jjj,17);         % time of ephemeris
    dt  = t_sec(ttt) - toe;          % seconds difference from epoch
    a   = sat_ephem(jjj,5)^2;        % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = sat_ephem(jjj,4);          % eccentricity
    n0  = sqrt(muE / a^3);           % nominal mean motion (rad/s)
    n   = n0 + sat_ephem(jjj,3);     % corrected mean motion, delta_n = gps_ephem_all(:,3)
    M   = sat_ephem(jjj,2) + n*dt;   % mean anomaly, M0 = gps_ephem_all(:,2)


    % Compute perigee, true and eccentric anomaly...
    % =====================================================================

    % Load argument of perigee to a local variable and add perigee rate, rad
    perigee  = sat_ephem(jjj,8) + sat_ephem(jjj,24) * dt;  

    % Compute Eccentric Anomaly, rad
    E = mean2eccentric(M,ecc);

    cosE = cos(E);  sinE = sin(E);

    % Compute rate of change of Eccentric Anomaly, rad/s (Eq. 20)
    E_dot = n / (1-ecc*cosE);

    % Compute true anomaly, rad
    nu = atan2(sqrt(1 - ecc*ecc).*sinE, ...
          (cosE - ecc)); 

    cosnu = cos(nu);  sinnu = sin(nu);  

    % Compute the argument of latitude, rad 
    phi = nu + perigee;  % true anomaly + argument of perigee


    % Compute corrections to argument of latitude, radius, and inclination
    % =====================================================================
    costwophi = cos(2*phi);  sintwophi = sin(2*phi);

    delta_u = sat_ephem(jjj,12) * sintwophi + ... % Cus = gps_ephem_all(jjj,12)
             sat_ephem(jjj,11) * costwophi;      % Cuc = gps_ephem_all(jjj,11)


    delta_r = sat_ephem(jjj,14) * sintwophi + ... % Crs = gps_ephem_all(jjj,14)
             sat_ephem(jjj,13) * costwophi;      % Crc = gps_ephem_all(jjj,13)

    delta_i = sat_ephem(jjj,16) * sintwophi + ... % Cis = gps_ephem_all(jjj,16)
             sat_ephem(jjj,15) * costwophi;      % Cic = gps_ephem_all(jjj,15)


    u = phi + delta_u;                                       % corrected argument of latitude
    r = a * (1 - ecc*cosE) + delta_r;                        % corrected radius  
    i = sat_ephem(jjj,7) + delta_i + sat_ephem(jjj,10) * dt; % corrected inclination 
                                                            % i_dot = sat_ephem(jjj,10)

    cosu = cos(u);  cos2u = cos(2*u);  
    sinu = sin(u);  sin2u = sin(2*u);


    % Compute Rates of Change of true anomaly, arg. of lat., radius, inclination
    % =====================================================================

    % Compute rate of change of true anomaly, rad/s  
    %   used in   http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
    nu_dot = sinE*E_dot*(1+ecc*cosnu) / (sinnu*(1-ecc*cosE));

    %Eq 24 and used in skyplot.cpp 
    %nu_dot2 = a*a*sqrt(1 - e^2)*n ./ (a .* (1 - e .* cos(E))).^2; 

    %the 2 previous equations for nu_dot are algebraically equal


    % Eq. 25, 26 and 24
    u_dot = nu_dot + 2*(sat_ephem(jjj,12)*cos2u-sat_ephem(jjj,11)*sin2u)*nu_dot;

    % Eq. 19, 20, 22 and 24 
    r_dot = a*ecc*sinE*n/(1-ecc*cosE) + 2*(sat_ephem(jjj,14)*cos2u-sat_ephem(jjj,13)*sin2u)*nu_dot;

    % Same format as Eq. 22 and 26 but with Cic and Cis instead
    i_dot = sat_ephem(jjj,10) + 2*(sat_ephem(jjj,16)*cos2u-sat_ephem(jjj,15)*sin2u)*nu_dot;

    i_dot2 = i_dot*i_dot; 



    % Compute satellite position in orbital plane (Eq. 13)
    % =====================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane



    % Compute satellite velocity in orbital plane, Eq. 18
    % =====================================================================
    xo_dot = r_dot*cosu - yo*u_dot;
    yo_dot = r_dot*sinu + xo*u_dot;



    % Corrected longitude of ascending node for node rate and Earth rotation
    % =====================================================================
    % Ascending node = ephem_all(jjj,6)
    node = sat_ephem(jjj,6) + (sat_ephem(jjj,9) - wE)*dt -  (wE * sat_ephem(jjj,17)); % Toe = gps_ephem_all(jjj,17)

    node_dot = sat_ephem(jjj,9) - wE;    %Eq. 10,  node rate = ephem_all(jjj,9)

    node_dot2 = node_dot*node_dot;


    % Calculate GPS Satellite Position in ECEF (m)
    % =====================================================================
    cosi = cos(i);     sini = sin(i);
    coso = cos(node);  sino = sin(node);


    % Satellite position in ECEF (m)
    x(ttt,1) = xo*coso - yo*cosi*sino;  %x-position  

    x(ttt,2) = xo*sino + yo*cosi*coso;  %y-position 

    x(ttt,3) = yo*sini;                 %z-position


    % Calculate Satellite Velocity in ECEF (m/s)
    % =====================================================================

    % Full velocity expression, Eq. 9
    %  Also presented in  http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
    v(ttt,1) = (xo_dot - yo*cosi*node_dot)*coso - (xo*node_dot + yo_dot*cosi - yo*sini*i_dot)*sino;

    v(ttt,2) = (xo_dot - yo*cosi*node_dot)*sino + (xo*node_dot + yo_dot*cosi - yo*sini*i_dot)*coso;

    v(ttt,3) = yo_dot*sini + yo*cosi*i_dot;



    % Keep track of health of each satellite
    % =====================================================================      
    health(ttt,1) = sat_ephem(jjj,25); % satellite health (0.00 is useable)





end % END of t_input loop =================================================
%==========================================================================    




    
    
