function [interp_sp3] = interp_sp3_orbits(sp3, dt, start_time, end_time)

% Re = 6378;

for prn = 1:32
	sp3_row_idx = sp3.data(:,sp3.col.PRN) == prn;

	sp3x = sp3.data(sp3_row_idx,sp3.col.X);
	sp3y = sp3.data(sp3_row_idx,sp3.col.Y);
	sp3z = sp3.data(sp3_row_idx,sp3.col.Z);
	sp3idx = sp3.data(sp3_row_idx,sp3.col.TOW);

	if isempty(sp3idx) % if no data exists for a prn
		continue
	end
	
	sp3_del_ts = sp3.data(sp3_row_idx,sp3.col.B)';
	
	tt = start_time:dt:end_time;
	
	if ( (~isempty(find(sp3x==0, 1))) || (~isempty(find(sp3x==0, 1)) ) || (~isempty(find(sp3x==0, 1))) )
		fprintf('sp3 zero val found');
	end
	
	for i = 1:length(tt)
		sp3pos = precise_orbit_interp(tt(i),sp3idx,...
			sp3x,...
			sp3y,...
			sp3z);
		
		interp_sp3.('x')(prn,i) = sp3pos(1);
		interp_sp3.('y')(prn,i) = sp3pos(2);
		interp_sp3.('z')(prn,i) = sp3pos(3);
		
		interp_sp3.('GPSSec')(prn,i) = tt(i);
	end
	%disp([num2str(prn) '|' num2str(length(sp3idx)) '|' num2str(length(tt)) '|']);
	
	interp_sp3.('clk_err')(prn,:) = interp1(sp3idx,sp3_del_ts,tt,'linear');
end

end

