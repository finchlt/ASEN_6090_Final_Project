function [WN2, TOW, WN1] = gpsvec2gpstow(gps_vec)

%==========================================================================
%==========================================================================
% [WN2, TOW, WN1] = gpsvec2gpstow(gps_vec)
%
% Calculates the GPS week number and seconds into the GPS week of a GPS
%  date/time vector.
%
%
% Author: Ben K. Bradley
% date: 04/28/2010
%
%
% INPUT:         Description                                      Units
%
%  gps_vec    - GPS date/time as a vector                     [y m d m h s]
%
%
% OUTPUT:       
%  
%  WN2        - week number count from 22Aug99 (rollover)          integer
%                This will be NaN if date is before 22Aug99
%  TOW        - seconds into the GPS week                          seconds
%  WN1        - week number count from 06jan80 (no rollover)       integer
%   
%
% Coupling:
%
%  Julian_Day
%
% References:
% 
%  [1] Montenbruck, O. and Gill, E. "Satellite Orbits: Models, Methods,
%          Applications". Corrected Third Printing, 2005. Springer-Verlag
%          Berlin Heidelberg. 2000.
%
%==========================================================================
%==========================================================================

yr=2000*(gps_vec(:,1)<20);
gps_vec(:,1)=gps_vec(:,1)+yr; % ONLY WORKS FOR year 2000+


% Store the fractional part of seconds
fracsec = gps_vec(6) - fix(gps_vec(6));
gps_vec(6)=fix(gps_vec(6));
    
% Set both GPS week number origins ========================================

JD_GPS_06jan80 = 2444244.5;  % Julian_Day([1980 1 6 0 0 0])

JD_GPS_22aug99 = 2451412.5;  % Julian_Day([1999 8 22 0 0 0])


% Error Check =============================================================
[JD_gps,JD1,JD2] = Julian_Day(gps_vec,0);


% Check if GPS date/time is before the start of GPS time (6jan1980 0hr)
if (JD_gps < 2444244.5)
    error('Input GPS date/time is before the start of GPS time (6jan1980 0hr)');
end


% Calculate Week Number ===================================================

WN1 = floor( (JD1+JD2 - JD_GPS_06jan80) / 7 );


if (JD_gps < JD_GPS_22aug99)
    
    WN2 = NaN;  
else

    WN2 = floor( (JD1+JD2 - JD_GPS_22aug99) / 7 );
end


% Calculate Time Of Week (TOW) ============================================
JD_Sunday = JD_GPS_06jan80 + 7*WN1;

TOW = round( (JD1+JD2 - JD_Sunday)*86400 ) + fracsec;  % seconds into GPS week








