function [range x_Tr t_Tt] = compute_range_sp3(sp3, SatPos, t, prn, userpos)

% Compute satellite postion
%[~,x,~] = broadcast2xv1(eph,t,prn);
x = SatPos;

sp3_row_idx = sp3.data(:,sp3.col.PRN) == prn;

sp3x = sp3.data(sp3_row_idx,sp3.col.X);
sp3y = sp3.data(sp3_row_idx,sp3.col.Y);
sp3z = sp3.data(sp3_row_idx,sp3.col.Z);
sp3idx = sp3.data(sp3_row_idx,sp3.col.TOW);

% Compute geomatric range
range = sqrt((x(1)-userpos(1))^2 + (x(2)-userpos(2))^2 + (x(3)-userpos(3))^2);

differ = 10000;

while (differ>1e-3)

range_old = range;

% compute the time of transmission
c = 299792458;
Tt = t(2) - range/c;
t_Tt = [1655 Tt];

% compute satellite position at Tt in ECEF
%[~,x_Tt,v,curreph] = broadcast2xv1(eph,t_Tt,prn);
x = precise_orbit_interp(Tt,sp3idx,...
			sp3x,...
			sp3y,...
			sp3z);

% rotate the satellite position to ECEF at Tr
w = 7.2921151467e-5;
phi = w * range/c;
RoMat = [cos(phi) sin(phi) 0;
         -sin(phi) cos(phi) 0;
         0 0 1];
x_Tr = RoMat * x*1e3;


% Compute new geometric range
range = sqrt((x_Tr(1)-userpos(1))^2 + (x_Tr(2)-userpos(2))^2 + (x_Tr(3)-userpos(3))^2);
differ = abs(range - range_old);

end


end
