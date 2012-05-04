%%% Generate/Plot SNR by track
%%% Praveen Vikram
tic

clear all
close all

station='p360';
year=2012;

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

truth = importdata([data_folder '/SnowDepth_p360.csv']);
snr_map = load([data_folder '/mat/snrMap']);

truth_arr = zeros(1,90);
snowIndex = zeros(1,90);
nTrack = zeros(1,90);
for doy = 1:90
	try
		filename=[station sprintf('%03d',doy) '0_' num2str(year) '_mp1.mat'];
		fprintf('=> loading Data File: %s | ',filename);
		clear data
		load([data_folder '/mat/' filename]);
	catch e
		disp(e)
		continue
	end
	
	% make truth array
	truth_idx = find(truth(:,1) == year & truth(:,2) == data.header.month & truth(:,3) == data.header.day);
	if ( (max(truth(truth(:,1) == year & truth(:,2) == data.header.month & truth(:,3) == data.header.day,12)) == 1) || ...
			(max(isnan(truth(truth(:,1) == year & truth(:,2) == data.header.month & truth(:,3) == data.header.day,12))) == 1) )
		truth_arr(doy) = 1;
		if isnan(truth(truth_idx,12))
			truth_arr(doy) = NaN;
		end
		if ( min(truth(truth_idx,12))==0 && max(truth(truth_idx,12)) == 1)
			truth_arr(doy) = 0.5;
		end
	end
	
	for prn = 1:32
		fprintf('%d,',prn);
		snr = data.S2(prn,:);
		elev = round(data.elev(prn,:));
		idx = find(elev>10);
		
		% check if snr val is within expected range
		for i = 1:length(idx)
			j=idx(i);
			if ~(  (snr(j)> (snr_map.mu(prn,elev(j))-3*snr_map.sig(prn,elev(j))^2) ) && ...
					 (snr(j)< (snr_map.mu(prn,elev(j))+3*snr_map.sig(prn,elev(j))^2) ) )
				 snowIndex(doy) = 1;
				 nTrack(doy) = nTrack(doy) + 1;
				 break;
			end
		end
	end
	disp('Done ')
	
end

sum(truth_arr == snowIndex)

return
%%
plot(1:90,truth_arr,'*')
hold on
plot(1:90,snowIndex,'ro')

xlabel('DOY')
legend('truth','Index from SNR')
title(['Snow Index for ' station ' - Year ' num2str(year)])
ylim([-0.5 1.5])
set(gca,'YTickLabel',['     ';'Clear';'Part ';'Snow '])