function [JD,JD1,JD2] = Julian_Day(tt,utcflag)

%==========================================================================
%==========================================================================
% [JD,JD1,JD2] = Julian_Day(time_vec)
%
% Calculates the Julian Day given a date/time vector. This algorithm is 
%  valid from the year 1900 to 2100 only. The Julian Day is output in its
%  total calculated form and as separate 0hr and fractional day components.
%  This is done for cases requiring high precision.  If utcflag is set to
%  one, then the leapsec.dat file is loaded to determine if the day is a
%  leap second day, in which case the Julian Day calculation uses a divisor 
%  of 61 instead of 60.
%
%
% Author: Ben K. Bradley
% date: 02/19/2009
%
%
% INPUT:           Description                                    Units
%
%  tt           - vector of date/time                         [y m d h m s]
%  utcflag      - flag stating if the input time is in UTC or not: 0=no,1=yes
% 
%
% OUTPUT:       
%    
%  JD           - Julian Day                     days since 1jan 12h 4713bc
%  JD1          - Julian Day of 0hr for input date
%  JD2          - hour, minute, second part of Julian Day
%                  e.g. JD = JD1 + JD2
%
%
% Coupling:
%
%   none
%
% References:
% 
%  [2] Vallado, D.A. "Fundamentals of Astrodynamics and Applications".
%         Third Edition. Microcosm Press. 2007.
%
%  [1] Meeus, J. "Astronomical Algorithms". Second Edition. With
%         Corrections as of June 15, 2005. Willmann-Bell, Inc. 1991 and
%         1998.
%
%==========================================================================
%==========================================================================


% Calculate Julian Day (1900 to 2100)

JD1 = 367*tt(1) - fix(1.75*(tt(1)+fix((tt(2)+9)/12))) + fix(275*tt(2)/9) + tt(3) + 1721013.5;


if (utcflag == 1)

    data = load('leapsec.dat');
    
    MJD0 = JD1 - 2400000.5;  % modified julian day of 0hr
    
    leap = 60 + sum( (MJD0==data(:,4)-1) );
    
else
    leap = 60;
end


JD2 = ((tt(6)/leap + tt(5))/60 + tt(4))/24;

JD = JD1 + JD2;




%Meeus p61 -----------------------------------------------------------

% if (Y > 1582) %Gregorian Calender
%     
%     A = fix(Y/100);
%     B = 2 - A + fix(A/4);
%     
% elseif ((Y+Mo/12+Day/31) <= (1582+10/12+4/31)) %Julian Calender
%     
%     B = 0;
%     
% else %Gregorian Calender
%     
%     A = fix(Y/100);
%     B = 2 - A + fix(A/4);
% end
% 
% 
% 
% if (Mo <= 2)
%     Y  = Y - 1;
%     Mo = Mo + 12;
% end
% 
% 
% 
% JD = fix(365.25*(Y + 4716)) + fix(30.6001*(Mo + 1)) ...
%      + Day + B - 1524.5 + ((Sec/60 + Mi) / 60 + Hr) / 24;
    
    
   
    
    
    