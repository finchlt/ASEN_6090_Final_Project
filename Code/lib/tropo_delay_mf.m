function [tropo_del] = tropo_delay_mf(n, tz, elev)
% elevation should be in deg

switch(n)
	case 1
		mf = 1/sind(elev);
		tropo_del = tz*mf;
	case 2
		mf = 1/sqrt(1 - (cosd(elev)/1.001)^2);
		tropo_del = tz*mf;
	otherwise
		tropo_del = [];
		disp('Invalid mapping function')
		return
end
		
	
end