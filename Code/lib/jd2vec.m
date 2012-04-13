function [time] = jd2vec(JD)

%==========================================================================
%==========================================================================
% [time] = jd2vec(JD)
%
% Converts a Julian date into a Gregorian calender date/time. This 
%  algorithm is accurate to 1e-4 seconds.
%
%
% Author: Ben K. Bradley
% date: 06/23/2009
%
%
% INPUT:        Description                                     Units
%
%  JD         - Julian date                   days since noon Jan. 1st 4713bc
%
% OUTPUT:       
%    
%   time      - date/time vector                              [y m d h m s]
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


% Calculate Modified Julian Date ==========================================
MJD = JD - 2400000.5;

% Julian Date at noon is calculated from MJD
a = floor(MJD) + 2400001;

% Fraction of day
q = MJD - floor(MJD);

% Auxiliary quantities, b,c,d,e,f =========================================
if (a < 2299161)
    b = 0;
else
    b = floor((a - 1867216.25)/36524.25);
end

if (a < 2299161)
    c = a + 1524;
else
    c = a + b - floor(b/4) + 1525;
end

d = floor((c - 121.1)/365.25);

e = floor(365.25*d);

f = floor((c - e)/30.6001);

% Calculate Day of Month ==================================================
D = c - e - floor(30.6001*f);

% Calculate Month of the year
Mo = f - 1 - 12*floor(f/14);

% Calculate Year
Y = d - 4715 - floor((7 + Mo)/10);


% Convert q into seconds into the day
q = q*86400;

% Calculate hour
H = fix(q / 3600);
q = q - H*3600;

% Calculate minute
Mi = fix(q / 60);
q = q - Mi*60;

% Calculate seconds
Sec = q;


% Assemble date/time vector ===============================================
time = [Y Mo D H Mi Sec];





