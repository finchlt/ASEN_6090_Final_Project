% ASEN 6090
% Daily MP RMS Reader

function []=plot_rms_mp(filename)
%% Load data file
data=load(filename);

%% Doy
doy=data(:,1)+(data(:,2)/366);
doy=data(:,2);
mp=data(:,3);
smoothmp=smooth(mp);
%% Create MP1 plot

figure
plot(doy,mp)
hold on
plot(doy,smoothmp,'r')
xlabel('DOY')
ylabel('MP1 [m]')
legend('RMS MP1','Smoothed RMS MP1')
figure
plot(doy,mp-smoothmp)
xlabel('DOY')
ylabel('MP1 Difference [m]')
