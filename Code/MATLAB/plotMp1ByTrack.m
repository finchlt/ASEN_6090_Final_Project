%%% Generate/Plot MP1
%%% Praveen Vikram
tic

clear all
close all

station='ab33';
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

prn = 17;

snow_fig=figure; hold on
clear_fig=figure; hold on
for doy = 33:90
	try
		filename=[station sprintf('%03d',doy) '0_' num2str(year) '_mp1.mat'];
		fprintf('=> loading Data File: %s | ',filename);
		clear data
		load([data_folder '/mat/' filename]);
		
		if(isfield(data.MP1,['PRN' num2str(prn)]))
			mp1 = data.MP1.(['PRN' num2str(prn)]);
		else
			continue
		end
		if max(truth(truth(:,1) == year & truth(:,2) == data.header.month & truth(:,3) == data.header.day,12)) == 1
			col = 'ro';
			disp('snow day | ')
			figure(snow_fig),plot(mp1(2,:),mp1(1,:),col,'MarkerSize',10)
		else
			col = 'go';
			figure(clear_fig),plot(mp1(2,:),mp1(1,:),col,'MarkerSize',10)
		end
		
		
		disp('Done ')
	catch e
		disp(e)
		continue
	end
end

figure(snow_fig),
title(['Antenna Covered By Snow [PRN ' num2str(prn) ']'])
xlabel('Elevation')
ylabel('MP1')
ylim([-5 5])
xlim([10 90])
grid on

figure(clear_fig),
title(['No Snow on Antenna [PRN ' num2str(prn) ']'])
xlabel('Elevation')
ylabel('MP1')
ylim([-5 5])
xlim([10 90])
grid on
return
%%

prn=2;
mp1 = data.MP1.(['PRN' num2str(prn)]);
figure,plot(mp1(2,:),mp1(1,:))
title(['PRN ' num2str(prn)])
xlabel('Elevation')
ylabel('MP1')
ylim([-5 5])
xlim([10 90])
grid on