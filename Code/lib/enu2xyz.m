function [XYZ] = enu2xyz(refLat, refLong, refH, e, n, u)
  % Convert east, north, up coordinates (labeled e, n, u) to ECEF
  % coordinates. The reference point (phi, lambda, h) must be given. 
  % All distances are in metres and angles in degrees
  %
  % Code taken from http://wiki.gis.com/wiki/index.php/Geodetic_system
 
  [XYZr] = wgslla2xyz(refLat,refLong, refH); % location of reference point
  Xr=XYZr(1);
  Yr=XYZr(2);
  Zr=XYZr(3);
  phiP = atan2(Zr,sqrt(Xr^2+Yr^2)); % Geocentric latitude
 
  X = -sin(refLong)*e - cos(refLong)*sin(phiP)*n + cos(refLong)*cos(phiP)*u + Xr;
  Y =  cos(refLong)*e - sin(refLong)*sin(phiP)*n + cos(phiP)*sin(refLong)*u + Yr;
  Z = cos(phiP)*n + sin(phiP)*u + Zr;
  
  XYZ=[X;Y;Z];