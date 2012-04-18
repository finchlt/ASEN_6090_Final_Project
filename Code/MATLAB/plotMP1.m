%%% Generate/Plot MP1
%%% Praveen Vikram

clear all
close all

station='p360';
doy='053';
year='2012';

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

filename=[station doy '0.mat'];

addpath(lib_folder);

data=load([data_folder '/' filename]);

% Check if P1 or C1 should be used
if ( (max(data.P1(:,1)) == 0) && (max(data.C1(:,1)) ~= 0) )
	rhoC=data.C1;
elseif (max(data.P1(:,1)) ~= 0)
	rhoC=data.P1;
end

% L1PhaseAdjust=zeros(1,32);
% L2PhaseAdjust=zeros(1,32);
% 
% L1 = zeros(32,size(rho,2));
% L2 = zeros(32,size(rho,2));
% 
% for prn = 1:32
% 	for epoch = 1:size(rho,2)
% 		if rho(prn,epoch) ~= 0
% 			if L1PhaseAdjust(prn) == 0
% 				L1PhaseAdjust(prn) = round((rho(prn,epoch)-lambda_L1*data.L1(prn,epoch))/lambda_L1);
% 				L2PhaseAdjust(prn) = round((rho(prn,epoch)-lambda_L2*data.L2(prn,epoch))/lambda_L2);
% 			end
% 			L1(prn,epoch) = data.L1(prn,epoch) + L1PhaseAdjust(prn);
% 			L2(prn,epoch) = data.L2(prn,epoch) + L2PhaseAdjust(prn);
% 		end
% 	end
% end

prn = 17;
si=1216;
ei=2822;

rho=rhoC(prn,si:ei);
L1=data.L1(prn,si:ei);
L2=data.L2(prn,si:ei);

L1PhaseAdjust = round((rho(1)-lambda_L1*L1(1))/lambda_L1);
L2PhaseAdjust = round((rho(1)-lambda_L2*L2(1))/lambda_L2);

L1 = L1 + L1PhaseAdjust;
L2 = L2 + L2PhaseAdjust;

mp1 = MP1(rho,L1,L2);

plot(mp1)

return
%%
for prn=1:32
plot(data.L1(prn,:)), hold on
end
