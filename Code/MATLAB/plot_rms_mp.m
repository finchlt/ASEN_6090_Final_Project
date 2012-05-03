% ASEN 6090
% Daily MP RMS Reader for ab33 p360 p455

%% Load data files
clear all;close all;clc;
dataall{1}=[];dataall{2}=[];dataall{3}=[];
doy{1}=[];doy{2}=[];doy{3}=[];

datafolder='../../Data/mp/';
stationname={'ab33','p360','p455'};
for s=1:3  % all 3 sites
    %% Load Data Files
    for k=9:12
        filename=sprintf('%s%s.mp.%02d',datafolder,stationname{s},k);
        data{s,k}=load(filename);
        dataall{s}=[dataall{s};data{s,k}];
        doy{s}=[doy{s};data{s,k}(:,1)+(data{s,k}(:,2)/365)];
    end
    
    %% Sort Data
    mp{s}=dataall{s}(:,3);
    smoothmp{s}=smooth(mp{s});
    diffmp{s}=mp{s}-smoothmp{s};
    stddiff{s}=std(diffmp{s});
    indmpout=abs(diffmp{s}-mean(diffmp{s}))>stddiff{s}*2;
    diffmpout_s2{s}=diffmp{s}(indmpout,:);
    doympout_s2{s}=doy{s}(indmpout,:);
    indmpout=abs(diffmp{s}-mean(diffmp{s}))>stddiff{s}*3;
    diffmpout_s3{s}=diffmp{s}(indmpout,:);
    doympout_s3{s}=doy{s}(indmpout,:);
    %% Create Plots
    figure
    plot(doy{s},mp{s})
    hold on
    plot(doy{s},smoothmp{s},'r')
    title(strcat(upper(stationname{s}),' RMS MP1'))
    xlabel(' Year')
    ylabel('RMS MP1 [m]')
    legend('RMS MP1','Smoothed')
    
    figure
    hold on
    plot(doy{s},diffmp{s})
    plot(doy{s},smooth(diffmp{s}),'k')
    plot(doympout_s2{s},diffmpout_s2{s},'ro')
    plot(doympout_s3{s},diffmpout_s3{s},'go')
    title(strcat(upper(stationname{s}),' RMS MP1 Difference from Smoothed'))
    xlabel(' Year')
    ylabel('RMS MP1 Difference')
    legend('RMS MP1 Difference','2 \sigma','3 \sigma')
    
end % End s=1:3  site for loop

%% Doy
% doy=data(:,1)+(data(:,2)/366);
% doy=data(:,2);
% mp=data(:,3);
% smoothmp=smooth(mp);
%% Create MP1 plot
% 
% figure
% plot(doy,mp)
% hold on
% plot(doy,smoothmp,'r')
% xlabel('DOY')
% ylabel('MP1 [m]')
% legend('RMS MP1','Smoothed RMS MP1')
% figure
% plot(doy,mp-smoothmp)
% xlabel('DOY')
% ylabel('MP1 Difference [m]')
