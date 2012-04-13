%%% Plot elevetation of all tracks
%%% Praveen Vikram

clear all
close all
% clc;

cur_folder = [pwd '/'];
lib_folder = [pwd '/../lib'];
data_folder = [pwd '/../../Data'];

addpath(lib_folder);

rinex_obs_file = 'p3600530.12o';
rinex_nav_file = 'auto0530.12n';

cutoff_elev = 10; % degrees

% some constants
c = 299792458; %speed of light in m/s
omegaE = 7.2921151467e-5; % 2*pi/86400; Earth rotation Rate
f_L1 = 1.57542e3;         % L1 freq in MHz
f_L2 = 1.22760e3;         % L2 freq in MHz
lambda_L1 = c/(f_L1*1e6); % L1 wavelength in m
lambda_L2 = c/(f_L2*1e6); % L2 wavelength in m

%% Read in Data
% 
% % Read Observation File
% fprintf('=> Using Rinex Observable File: %s',rinex_obs_file);
% 
% % Obtain site xyz from Rinex header 
% % and observables
% [rinex_obs_fid, site_xyz, data, ~] = read_rinex_obs([data_folder '/' rinex_obs_file]);
% 
% if rinex_obs_fid == -1
% 	fprintf('... Not Found!\n');
% 	fprintf('\n ==> Aborting!! <==\n');
% 	return
% else
% 	fprintf('\t=>Done!\n');
% end
% 
% % Read Nav File
% fprintf('=> Using Rinex Nav File: %s',rinex_nav_file);
% [gps_ephem,ionoparams] = read_GPSbroadcast([data_folder '/' rinex_nav_file]);
% fprintf('... Done!\n');
load p3600520.mat

%% Fun Part

% transfer user position from ECEF to LLA
[wlat, wlon, walt] = wgsxyz2lla(site_xyz);

N_epochs = length(data.GPSSec);

% Generate some vars to store data
data.('elev') = zeros(32,1);
data.('azi') = zeros(32,1);

for epoch = 1:N_epochs
	
	visible_prn = find(data.P2(:,epoch) > 0);
	N_obs = length(visible_prn);
	
	%fprintf(' %03d,',epoch);
	fprintf('.');
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
		
		%compute elevation angle
		elev = 90-acos(enu(3,1)/range_enu)*180/pi;

		data.('elev')(prn,epoch) = elev;

	end
	
end

data.elev(data.elev == 0) = NaN;
data.azi(data.azi == 0) = NaN;

disp('Done')

return
%%

figure, hold on

l=[];
c=['r-o';'g-o';'b-o';'c-o';'m-o';'y-o';'k-o';...
	'r-*';'g-*';'b-*';'c-*';'m-*';'y-*';'k-*';...
	'r--';'g--';'b--';'c--';'m--';'y--';'k--';...
	'r-.';'g-.';'b-.';'c-.';'m-.';'y-.';'k-.';];

for prn=1:32
	disp(num2str(prn))
	a=result.azi(prn,:);
	e=result.elev(prn,:);

	plot(e,c(prn,:));
	waitforbuttonpress

end

%%

prn = 17
disp(num2str(prn))
a=result.azi(prn,:);
e=result.elev(prn,:);

plot(e);

