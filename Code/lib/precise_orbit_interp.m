function [XYZpos] = precise_orbit_interp(Ttime, epochs, X, Y, Z)

sidereal_day = 0.99726956634;
period = sidereal_day;

P0 = 2.0*pi / period;

Time = epochs;
% keyboard
%find the nearest zero point
[~, closestI] = min(abs(epochs - Ttime));

NDAT = 4;
NDAT1 = NDAT+1;
I2NDAT = 2*NDAT;
I2NDAT1 = I2NDAT + 1;

% data popints to use
if closestI < NDAT1
    k=[1:I2NDAT1];
end
if closestI>length(epochs)-NDAT
    k=[length(epochs)-I2NDAT:length(epochs)];
end
if closestI>=NDAT1&&closestI<=length(epochs)-NDAT
    k = [closestI-NDAT:closestI+NDAT];
end

Xi = X(k);
Yi = Y(k);
Zi = Z(k);
Timei = (Time(k)-Ttime)/86400;

Nest = I2NDAT1;
Nmeas = length(k);
A  = zeros(Nmeas, Nest);
A(:,1) = ones(Nmeas,1); 
B  = zeros(1, Nest);
B(1,1) = 1;

gps_rel_time = 0;

for i = 1:NDAT
  kk = 2+(i-1)*2;
  P =  P0*i;
  A(:,kk) = sin(P*Timei);
  A(:,kk+1) = cos(P*Timei);
  B(1,kk) = sin(P*gps_rel_time);
  B(1,kk+1) = cos(P*gps_rel_time);
end

XCoeffs = A\Xi; 
YCoeffs = A\Yi; 
ZCoeffs = A\Zi; 

XYZpos = [B*XCoeffs B*YCoeffs B*ZCoeffs]';
% XYZpos = [newX; newY; newZ]*1e6;
%keyboard
end
