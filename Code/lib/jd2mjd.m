function [MJD] = jd2mjd(JD)

%==========================================================================
%==========================================================================
% [MJD] = jd2mjd(JD)
%
% Calculates the Modified Julian Day from a Julian Day 
%
%
% Author: Ben K. Bradley
% date: 06/24/2009
%
%
% INPUT:         Description                                         Units
%
%  JD       - Julian Day                         days since 1jan 12h 4713bc
%
% OUTPUT:       
%    
%  MJD      - Modified Julian Day                 days since 17nov 0h 1858
%
%
% Coupling:
%
%   none
%
% References:
% 
%   [1] Montenbruck, O. and Gill, E. "Satellite Orbits: Models, Methods,
%          Applications". Corrected Third Printing, 2005. Springer-Verlag
%          Berlin Heidelberg. 2000.
%
%==========================================================================
%==========================================================================


% Calculate Modified Julian Day 
MJD = JD - 2400000.5;



