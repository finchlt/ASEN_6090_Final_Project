% ASEN 6090
% Daily MP RMS Reader for ab33 p360 p455

%% Load data files
clear all;close all;clc;
dataall{1}=[];dataall{2}=[];dataall{3}=[];
snowdataall{1}=[];snowdataall{2}=[];snowdataall{3}=[];
doy{1}=[];doy{2}=[];doy{3}=[];
doysnow{1}=[];doysnow{2}=[];doysnow{3}=[];

datafolder='../../Data/';
stationname={'ab33','p360','p455'};
for s=1:3  % all 3 sites
    %% Load Data Files
    for k=9:12
        filename=sprintf('%smp/%s.mp.%02d',datafolder,stationname{s},k);
        data{s,k}=load(filename);
        dataall{s}=[dataall{s};data{s,k}];
        doy{s}=[doy{s};data{s,k}(:,1)+(data{s,k}(:,2)/365)];
        snowfile=sprintf('%ssnow/%s_20%02d.csv',datafolder,stationname{s},k);
        snowdata{s,k}=load(snowfile);
        snowdataall{s}=[snowdataall{s};snowdata{s,k}];
        doysnow{s}=[doysnow{s};snowdata{s,k}(:,2)+(snowdata{s,k}(:,3)/365)];
    end
    
    %% Sort Data
    mp{s}=dataall{s}(:,3);
    smoothmp{s}=smooth(mp{s});
    diffmp{s}=mp{s}-smoothmp{s};
    stddiff{s}=std(diffmp{s});
    
    % 2*std
    ind2sig=abs(diffmp{s}-mean(diffmp{s}))>stddiff{s}*2;
    diffmp_s2{s}=diffmp{s}(ind2sig,:);
    doymp_s2{s}=doy{s}(ind2sig,:);
    doysnow_s2{s}=doysnow{s}(ind2sig,:);
    snow_s2{s}=snowdataall{s}(ind2sig,5);
    
    % 3*std
    ind3sig=abs(diffmp{s}-mean(diffmp{s}))>stddiff{s}*3;
    diffmp_s3{s}=diffmp{s}(ind3sig,:);
    doymp_s3{s}=doy{s}(ind3sig,:);
    doysnow_s3{s}=doysnow{s}(ind3sig,:);
    snow_s3{s}=snowdataall{s}(ind3sig,5);
    
    %% Create Plots
    figure
    plot(doy{s},mp{s},'linewidth',1.25)
    hold on
    grid on
    plot(doy{s},smoothmp{s},'r','linewidth',1.25)
    title(strcat(upper(stationname{s}),' RMS MP1'),'fontsize',16)
    xlabel('Year','fontsize',14)
    ylabel('RMS MP1 [m]','fontsize',14)
    legend('RMS MP1','Smoothed')
    % Print
    h=gcf;
    str=sprintf('Figures/%s_mp',stationname{s});
    print(h,'-dpng',str);
    
    % MP1 Difference
    figure
    hold on
    grid on
    plot(doy{s},diffmp{s},'linewidth',1.25)
    plot(doymp_s2{s},diffmp_s2{s},'ro','markersize',4,'linewidth',3)
    plot(doymp_s3{s},diffmp_s3{s},'ko','markersize',4,'linewidth',3)
    title(strcat(upper(stationname{s}),' RMS MP1 Difference from Smoothed'),'fontsize',16)
    xlabel('Year','fontsize',14)
    ylabel('RMS MP1 Difference [m]','fontsize',14)
    legend('RMS MP1 Difference','2 \sigma','3 \sigma')
    % Print
    h=gcf;
    str=sprintf('Figures/%s_mpdiff',stationname{s});
    print(h,'-dpng',str);
    
    
    figure
    hold on
    grid on
    plot(doysnow{s},snowdataall{s}(:,5),'linewidth',1.25)
    plot(doysnow_s2{s},snow_s2{s},'ro','markersize',4,'linewidth',3)
    plot(doysnow_s3{s},snow_s3{s},'ko','markersize',4,'linewidth',3)
    xlabel('Year','fontsize',14)
    ylabel('Snow Depth [m]','fontsize',14)
    legend('Snow Depth','3 \sigma Outliers')
    title(strcat(upper(stationname{s}),' Outliers and Snow Depth'),'fontsize',16)
    % Print
    h=gcf;
    str=sprintf('Figures/%s_snow',stationname{s});
    print(h,'-dpng',str);
    
    
end % End s=1:3  site for loop

%% Compare P360 to Truth
% Load Truth Data
% p360truth=load()

% Trim p360 data to 2012+-120
trimdoy=(doy{2}-2012)*365;
inddoy=(trimdoy>-120);
trimdoy=trimdoy(inddoy);
trimdiffmp=diffmp{2}(inddoy);
% trim truthdiffmp and truthdoy
% truthdoy=trimdoy(logical(truth),:)
% truthdiffmp=trimdiffmp(logical(truth),:)

% Trim 2 and 3 sigma data
trimdoy_s2=(doysnow_s2{2}-2012)*365;
inds2=trimdoy_s2>-120;
trimdoy_s2=trimdoy_s2(inds2);
trimdoy_s3=(doysnow_s3{2}-2012)*365;
inds3=trimdoy_s3>-120;
trimdoy_s3=trimdoy_s3(inds3);

% Plot MP
figure
hold on
plot(trimdoy,trimdiffmp,'linewidth',1.25)
% plot(truthdoy,truthdiffmp,'ko','markersize',4,'linewidth',3)
xlabel('Days Since 2012','fontsize',14)
ylabel('RMS MP1 Difference [m]','fontsize',14)
legend('MP1 Difference','Snow Days')

% plot 2 and 3 sigma data
figure
hold on
% plot(p360truth,p360truth)
plot(trimdoy_s2,.9,'ro','markersize',4,'linewidth',3)
plot(trimdoy_s3,.95,'ko','markersize',4,'linewidth',3)
ylabel('Flag','fontsize',14)
xlabel('Days Since 2012','fontsize',14)
title('P360 Truth Flags and Outliers','fontsize',16)
legend('P360 Truth','2 \sigma','3 \sigma')
h=gcf;
str=sprintf('Figures/p360_truth');
print(h,'-dpng',str);