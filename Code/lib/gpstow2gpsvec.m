function [GPS_vec] = gpstow2gpsvec(WN,TOW,rollflag)

%==========================================================================
%==========================================================================
% [GPS_vec] = gpstow2gpsvec(WN,TOW,rollflag)
%
% Calculates a GPS date/time vector from the week number and time of week.
%
%
% Author: Ben K. Bradley
% date: 04/28/2010
%
%
% INPUT:         Description                                        Units
%
%  WN         - GPS week number                                    integer
%  TOW        - seconds into the GPS week                          seconds
%                (starting from 0hr on Sunday)
%  rollflag   - dictates week number origin of input            1 = 06jan80
%                                                               2 = 22Aug99
%
% OUTPUT:       
%  
%  gps_vec    - GPS date/time as a vector                     [y m d m h s]
%
%
% Coupling:
%
%  jd2vec
%
% References:
% 
%  [1] Montenbruck, O. and Gill, E. "Satellite Orbits: Models, Methods,
%          Applications". Corrected Third Printing, 2005. Springer-Verlag
%          Berlin Heidelberg. 2000.
%
%==========================================================================
%==========================================================================

fracsec = TOW - fix(TOW);


% Set both GPS week number origins ========================================

%JD_GPS_06jan80 = 2444244.5;  % Julian_Day([1980 1 6 0 0 0])

%JD_GPS_22aug99 = 2451412.5;  % Julian_Day([1999 8 22 0 0 0])

JD_GPS_epoch = [2444244.5   2451412.5];


% Compute Julian Date  ====================================================
    
JD_gps = WN*7 + JD_GPS_epoch(rollflag) + (TOW/86400);

    
% Convert Julian Date to date/time vector (Gregorian) =====================
    
GPS_vec = jd2vec(JD_gps);


GPS_vec(6) = fix(GPS_vec(6)) + fracsec;















