%% SNR Plots
% Produce Plots of the SNR data for each of the two different days
% Logan Finch
clc; clear all; close all;

%% Load Data
load('snr_data.mat');

%% Sort by PRN
prn = 17;

snr21_index = p360_feb21.col.S1;
snr22_index = p360_feb21.col.S1;

index21 = find(p360_feb21.data(:,3) == prn);
index22 = find(p360_feb22.data(:,3) == prn);

p360_feb21.data = p360_feb21.data(index21,:);
p360_feb22.data = p360_feb22.data(index22,:);

S1_21 = p360_feb21.data(:,snr21_index);
S1_22 = p360_feb22.data(:,snr22_index);

epoch21 = p360_feb21.data(:,2);
epoch22 = p360_feb22.data(:,2);

lla = ecef2lla(rec_pos_p360');


for i = 1 : length(epoch21)
    x_21 = SatPos(sp3_21, epoch21(i), 1,prn);
    [~,el21(i,1),~] = ecef2azelrange2(x_21,rec_pos_p360,lla(1),lla(2));
end

for i = 1:length(epoch22);
    x_22 = SatPos(sp3_22, epoch22(i), 1,prn);
    [~,el22(i,1),~] = ecef2azelrange2(x_22,rec_pos_p360,lla(1),lla(2));
end

%% Plot Data
sidereal = (24-23.9344696);
tp21 = (epoch21-epoch21(1))./3600;
tp22 = (epoch22-epoch22(1))./3600-sidereal;
tp22 = tp22-tp22(1);

ind21 = find(tp21<7);
ind22 = find(tp22<7);

el21sort = el21(ind21);
el22sort = el22(ind22);
tp21sort = tp21(ind21);
tp22sort = tp22(ind22);
S21_sort = S1_21(ind21,1);
S22_sort = S1_22(ind22,1);


figure
plot(tp21(ind21),smooth(S1_21(ind21,1)),'r')
hold on
plot(tp22(ind22),smooth(S1_22(ind22,1)),'b')
xlabel('Time (hours)','fontsize',14)
ylabel('SNR (dB-Hz)','fontsize',14)
title('SNR vs Time','fontsize',14)
legend('SNR Data for 2/21/12','SNR Data for 2/22/12')

figure
plot(el21sort,S1_21(ind21,1),'r')
hold on
plot(el22sort,S1_22(ind22,1),'b')
xlabel('Elivation Angle (degrees)','fontsize',14)
ylabel('SNR (dB-Hz)','fontsize',14)
title('SNR vs Elivation Angle','fontsize',14)
legend('SNR Data for 2/21/12','SNR Data for 2/22/12','location','southeast')

ind21 = find(el21sort>80);
ind22 = find(el21sort>80);

tp21high = tp21sort(ind21,1);
tp22high = tp22sort(ind22,1);
el21high = el21sort(ind21,1);
el22high = el22sort(ind22,1);
S21high = S21_sort(ind21,1);
S22high = S22_sort(ind22,1);


ind = 63;
figure
plot(el21high(1:ind),S21high(1:ind),'r')
hold on
plot(el22high(1:ind),S22high(1:ind),'b')
xlabel('Elivation Angle (degrees)','fontsize',14)
ylabel('SNR (dB-Hz)','fontsize',14)
title('SNR vs Elivation Angle','fontsize',14)
legend('SNR Data for 2/21/12','SNR Data for 2/22/12','location','southeast')


ind = 763;
figure
plot(el21sort(1:(ind-1)),smooth(S21_sort(1:(ind-1)),10)-smooth(S22_sort(2:ind),10))
