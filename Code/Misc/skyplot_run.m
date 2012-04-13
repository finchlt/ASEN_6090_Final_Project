%%% Make a SkyPlot
%%% Praveen Vikram

%%% FIX: Something messed up with Azi 

clear all
close all
% clc;

cur_folder = [pwd '/'];
lib_folder = [pwd '/../lib'];
data_folder = [pwd '/../../Data'];

addpath(lib_folder);

rinex_obs_file = 'p3600520.12o';
rinex_nav_file = 'auto0520.12n';

cutoff_elev = 10; % degrees

% some constants
c = 299792458; %speed of light in m/s
omegaE = 7.2921151467e-5; % 2*pi/86400; Earth rotation Rate
f_L1 = 1.57542e3;         % L1 freq in MHz
f_L2 = 1.22760e3;         % L2 freq in MHz
lambda_L1 = c/(f_L1*1e6); % L1 wavelength in m
lambda_L2 = c/(f_L2*1e6); % L2 wavelength in m

%% Read in Data

% Read Observation File
fprintf('=> Using Rinex Observable File: %s',rinex_obs_file);

% Obtain site xyz from Rinex header 
% and observables
[rinex_obs_fid, site_xyz, data, ~] = read_rinex_obs([data_folder '/' rinex_obs_file]);

if rinex_obs_fid == -1
	fprintf('... Not Found!\n');
	fprintf('\n ==> Aborting!! <==\n');
	return
else
	fprintf('\t=>Done!\n');
end

% Read Nav File
fprintf('=> Using Rinex Nav File: %s',rinex_nav_file);
[gps_ephem,ionoparams] = read_GPSbroadcast([data_folder '/' rinex_nav_file]);
fprintf('... Done!\n');

%% Fun Part

% transfer user position from ECEF to LLA
[wlat, wlon, walt] = wgsxyz2lla(site_xyz);

N_epochs = length(data.GPSSec);

% Generate some vars to store data
result.('elev') = zeros(32,1);
result.('azi') = zeros(32,1);

for epoch = 1:N_epochs
	
	visible_prn = find(data.P2(:,epoch) > 0);
	N_obs = length(visible_prn);
	
	fprintf(' %03d,',epoch);
	% Define variables
	
	geo_range = zeros(N_obs,1);
	SatPos = zeros(N_obs,3);
	elev = zeros(N_obs,1);
	azi = zeros(N_obs,1);
	
	data.N_obs(epoch) = N_obs;
	
	% Calculate the Iono-free combination (P3)
	for prn_idx = 1:N_obs
		prn = visible_prn(prn_idx);
		
		t = [data.GPSWeek(epoch) data.GPSSec(epoch)];
		
		%[geo_range(prn_idx) SatPos(prn_idx,:) t_Tt] = compute_range_sp3(sp3, SatPos, t, prn, site_xyz);
		[~, SatPos(prn_idx,:), ~, ~, curreph] = compute_range(gps_ephem, t, prn, site_xyz);
		
		%transfer from ECEF to ENU
		enu = wgsxyz2enu(SatPos(prn_idx,:)', wlat, wlon, walt);
		
		% compute range under enu
		range_enu = sqrt(enu'*enu);
				
		% compute elevation angle
		azi = atan2(SatPos(prn_idx,1),SatPos(prn_idx,2))*180/pi;
		if (azi<0) azi = azi+360; end
		
		%compute elevation angle
		elev = 90-acos(enu(3,1)/range_enu)*180/pi;

		result.('elev')(prn,epoch) = elev;
		result.('azi')(prn,epoch) = azi;
	end
	
end

result.elev(result.elev == 0) = NaN;
result.azi(result.azi == 0) = NaN;

return
%%

figure, hold on

for prn=2:2
	
	a=result.azi(prn,:);
	e=result.elev(prn,:);
% 	e=1-e;

	if exist('h') == 0
		h = skyplot(e,a,num2str(prn));
	else
		h = skyplot(e,a,num2str(prn),h);
	end
end
