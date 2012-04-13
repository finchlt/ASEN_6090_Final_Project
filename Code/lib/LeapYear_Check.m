function [leap] = LeapYear_Check(Yr)

%==========================================================================
%==========================================================================
% [leap] = LeapYear_Check(Yr) 
%
% Determines if given year is a leap year (i.e. Feb. has 29 days)
%  Leap years should be:  1980 1984 1988 1992 1996 2000 2004 2008 2012 2016
%
%
% Author: Ben K. Bradley
% date: 03/13/2009
%
%
% INPUT:          Description                                        Units
%
%  Yr         - year to be checked                                   years
% 
%
% OUTPUT:       
%    
%  leap       - true or false, Yr is a leap year          0 = not leap year
%                                                         1 = leap year
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
if (Yr <= 0)
    errordlg('Error in LeapYear_Check.m  :  Invalid year input');
end


%Check if current year is leap year
if ((rem(Yr,4)==0) && (rem(Yr,100)~=0)) || (rem(Yr,400)==0)
    %Year is leap year (Feb. has 29 days)
    leap=1;
else
    leap=0;
end




