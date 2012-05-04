%%% Generate/Plot MP1
%%% Praveen Vikram
tic

clear all
close all

station='p360';
year=2011;

% some constants
c = 299792458; %speed of light in m/s
omegaE = 7.2921151467e-5; % 2*pi/86400; Earth rotation Rate
f_L1 = 1.57542e3;         % L1 freq in MHz
f_L2 = 1.22760e3;         % L2 freq in MHz
lambda_L1 = c/(f_L1*1e6); % L1 wavelength in m
lambda_L2 = c/(f_L2*1e6); % L2 wavelength in m

cur_folder = [pwd '/'];
lib_folder = [pwd '/../lib'];
data_folder = [pwd '/../../Data'];

addpath(lib_folder);

for doy=200:249
	disp(['DOY: ' num2str(doy)])
	try
	
	% ================================================
	% calculate elevation angle for the entire dataset
	% ================================================
	yr_str=num2str(year);
	rinex_nav_file = ['auto' sprintf('%03d',doy) '0.' yr_str(3:4) 'n'];
	% Read Nav File
	fprintf('=> Using Rinex Nav File: %s',rinex_nav_file);
	[gps_ephem,ionoparams] = read_GPSbroadcast([data_folder '/rinex/' rinex_nav_file]);
	fprintf('... Done!\n');
	
	filename=[station sprintf('%03d',doy) '0_' num2str(year) '.mat'];
	fprintf('=> Using Data File: %s',filename);	
	data=load([data_folder '/mat/' filename]);
	fprintf('... Done!\n');
		
	N_epochs = length(data.sec);
	
	% Generate some vars to store data
	data.('elev') = zeros(32,1);
	
	for epoch = 1:N_epochs
		
		visible_prn = find(data.P2(:,epoch) > 0);
		N_obs = length(visible_prn);
		
		% fprintf(' %03d,',epoch);
		% fprintf('.');
		% Define variables
		
		geo_range = zeros(N_obs,1);
		SatPos = zeros(N_obs,3);
		elev = zeros(N_obs,1);
		
		data.N_obs(epoch) = N_obs;
		site_xyz = data.header.recvXYZ;
		[wlat, wlon, walt] = wgsxyz2lla(site_xyz);
		
		% Calculate the Iono-free combination (P3)
		for prn_idx = 1:N_obs
			prn = visible_prn(prn_idx);
			
			%t = [data.GPSWeek(epoch) data.GPSSec(epoch)];
			[~, month, day] = ydoy2ymd(year, doy);
			
			%cur_t = [data.header.year; data.header.month; data.header.day; data.hour(epoch); data.min(epoch); data.sec(epoch)]'; 
			%[~, gpssec, gpswk] = gpsvec2gpstow(cur_t);

			[gpswk gpssec] = GPSweek(year, month, day, data.hour(epoch), data.min(epoch), data.sec(epoch));
			t=[gpswk gpssec];
			%[geo_range(prn_idx) SatPos(prn_idx,:) t_Tt] = compute_range_sp3(sp3, SatPos, t, prn, site_xyz);
			[~, SatPos(prn_idx,:), ~, ~, curreph] = compute_range(gps_ephem, t, prn, site_xyz);
			
			%transfer from ECEF to ENU
			enu = wgsxyz2enu(SatPos(prn_idx,:)', wlat, wlon, walt);
			
			% compute range under enu
			range_enu = sqrt(enu'*enu);
			
			%compute elevation angle
			elev = 90-acos(enu(3,1)/range_enu)*180/pi;
			
			data.('elev')(prn,epoch) = elev;
			
		end
		if mod((epoch/N_epochs)*10,1) == 0 % print . every 10% of obs
			fprintf('.');
		end
	end
	disp('Done!')
	
	% =========================================
	% cycle through each prn and extract tracks
	% =========================================
	for prn = 1:32
		idx = find(data.elev(prn,:)>10);
		elev = data.elev(prn,idx);
		
		% check if P1 or C1 should be used
		if max(data.P1(prn,idx) ~= 0)
			rho = (data.P1(prn,:));
		elseif max(data.C1(prn,idx) ~= 0)
			rho = (data.C1(prn,:));
		else
			continue;
		end
		
		L1=data.L1(prn,:);
		L2=data.L2(prn,:);
		
		% need to phase adjust each seperate track 
		peaks = [0 find(diff(idx)>1) length(idx)];
		
		for i = 2:length(peaks)
		L1PhaseAdjust = round((rho(idx(peaks(i-1)+1))-lambda_L1*L1(idx(peaks(i-1)+1)))/lambda_L1);
		L2PhaseAdjust = round((rho(idx(peaks(i-1)+1))-lambda_L2*L2(idx(peaks(i-1)+1)))/lambda_L2);
		
		L1(idx(peaks(i-1)+1):idx(peaks(i))) = L1(idx(peaks(i-1)+1):idx(peaks(i))) + L1PhaseAdjust;
		L2(idx(peaks(i-1)+1):idx(peaks(i))) = L2(idx(peaks(i-1)+1):idx(peaks(i))) + L2PhaseAdjust;
		end
		mp1 = MP1(rho(idx),L1(idx),L2(idx));

		if(max(abs(mp1)) > 50)
			continue
		end
		
		data.('MP1').(['PRN' num2str(prn)]) = [mp1;elev]; 
	end
	
	% save the mp1 mat
	filename=[station sprintf('%03d',doy) '0_' num2str(year) '_mp1.mat'];
	fprintf('=> saving Data File: %s',filename);	
	save([data_folder '/mat/' filename],'data');
	disp(['Done! (' num2str(toc/60) 'mins)']);
	disp('-------------');
	catch e
		disp(e)
		continue
	end
end % end doy loop

return
%%

% Check if P1 or C1 should be used
if ( (max(data.P1(:,1)) == 0) && (max(data.C1(:,1)) ~= 0) )
	rhoC=data.C1;
elseif (max(data.P1(:,1)) ~= 0)
	rhoC=data.P1;
end

prn=17;shr=6;ehr=11;si=find(data.hour==shr,1,'first');ei=find(data.hour==ehr,1,'first');
%prn=14;si=4300;ei=5705;

rho=rhoC(prn,si:ei);
L1=data.L1(prn,si:ei);
L2=data.L2(prn,si:ei);

L1PhaseAdjust = round((rho(1)-lambda_L1*L1(1))/lambda_L1);
L2PhaseAdjust = round((rho(1)-lambda_L2*L2(1))/lambda_L2);

L1 = L1 + L1PhaseAdjust;
L2 = L2 + L2PhaseAdjust;

mp1 = MP1(rho,L1,L2);

%=======

doy='052';
year='2012';

filename=[station doy '0.mat'];

addpath(lib_folder);

data=load([data_folder '/' filename]);

% Check if P1 or C1 should be used
if ( (max(data.P1(:,1)) == 0) && (max(data.C1(:,1)) ~= 0) )
	rhoC=data.C1;
elseif (max(data.P1(:,1)) ~= 0)
	rhoC=data.P1;
end

prn=17;shr=6;ehr=11;si=find(data.hour==shr,1,'first')+4;ei=find(data.hour==ehr,1,'first')+4;
%prn=14;si=4300;ei=5705;

rho=rhoC(prn,si:ei);
L1=data.L1(prn,si:ei);
L2=data.L2(prn,si:ei);

L1PhaseAdjust = round((rho(1)-lambda_L1*L1(1))/lambda_L1);
L2PhaseAdjust = round((rho(1)-lambda_L2*L2(1))/lambda_L2);

L1 = L1 + L1PhaseAdjust;
L2 = L2 + L2PhaseAdjust;

mp1_2 = MP1(rho,L1,L2);

%%

plot(smooth(mp1,10))
hold on
plot(smooth(mp1_2,10),'r')
legend('doy 053','doy 052')
ylabel('MP_1')
xlabel('epoch')
title([station ' - ' num2str(prn) ' - ' num2str(doy) '/' num2str(year)])
% ylim([0 1])
return
%%
for prn=1:32
plot(data.L1(prn,:)), hold on
end
