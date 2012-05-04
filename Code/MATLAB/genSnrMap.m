%%% Generate SNR map
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

%truth = importdata([data_folder '/SnowDepth_p360.csv']);
mu = zeros(32,90);
interim = zeros(32,90);
n = zeros(32,90);
sig = zeros(32,90);

for doy = 200:250
	filename=[station sprintf('%03d',doy) '0_' num2str(year) '_mp1.mat'];
	fprintf('=> loading Data File: %s | ',filename);
	try
		clear data
		load([data_folder '/mat/' filename]);
	catch
		disp('Failed')
		continue
	end
	for prn = 1:32
		snr = data.S2(prn,:);
		elev = round(data.elev(prn,:));
		idx = find(elev>10);
		
		for i=1:length(idx)
			j = idx(i);
			elev_idx=elev(j);
			
			mu(prn,elev_idx) = (n(prn,elev_idx)*mu(prn,elev_idx) + snr(j))/(n(prn,elev_idx)+1);
			interim(prn,elev_idx) = (n(prn,elev_idx)*interim(prn,elev_idx) + (snr(j))^2)/(n(prn,elev_idx)+1);
			n(prn,elev_idx) = n(prn,elev_idx) + 1;
			
			sig(prn,elev_idx) = sqrt(interim(prn,elev_idx)-(mu(prn,elev_idx))^2);
		end
		fprintf('%d,',prn)
	end
	disp('Done ')
end

save([data_folder '/mat/snrMap'], 'mu', 'sig', 'n') 
return

%%
col = ['r.';'g.';'b.';'c.';'m.';'y.';...
	'rd';'gd';'bd';'cd';'md';'yd';...
	'ro';'go';'bo';'co';'mo';'yo';...
	'r*';'g*';'b*';'c*';'m*';'y*';...
	'rs';'gs';'bs';'cs';'ms';'ys';...
	'rp';'gp';'bp';'cp';'mp';'yp'];
for prn = 1:32
errorbar(1:90,mu(prn,:),sig(prn,:),col(prn,:))
hold on
end

ylim([15 55])
xlim([10 90])
xlabel('Elevation')
ylabel('SNR')
title(['SNR Map for ' station])
