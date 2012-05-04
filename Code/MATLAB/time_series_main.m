%% Time Series Data Analysis

clc; clear all; close all;

%% AB33
sd = [2009,1,1];
ed = [2012,4,1];
figure
subplot(2,1,1)
TimeSeriesPlot('AB33_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('ab33',sd,ed)
str = sprintf('AB33_long.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..

sd = [2011,9,3];
ed = [2012,4,29];
figure
subplot(2,1,1)
TimeSeriesPlot('AB33_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('ab33',sd,ed)
str = sprintf('AB33_2012.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..
TimeSeriesYearPlot('AB33_timeSeries.mat')
h = gcf;
cd Figures/
print(h,'-dpng','AB33_doy.png')
cd ..


%% P360
sd = [2009,1,1];
ed = [2012,4,1];
figure
subplot(2,1,1)
TimeSeriesPlot('P360_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('p360',sd,ed)
str = sprintf('P360_long.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..

sd = [2011,9,3];
ed = [2012,4,29];
figure
subplot(2,1,1)
TimeSeriesPlot('P360_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('p360',sd,ed)
str = sprintf('P360_2012.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..

TimeSeriesYearPlot('P360_timeSeries.mat')
h = gcf;
cd Figures/
print(h,'-dpng','P360_doy.png')
cd ..

%% P455
sd = [2009,1,1];
ed = [2012,4,1];
figure
subplot(2,1,1)
TimeSeriesPlot('P455_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('p455',sd,ed)
str = sprintf('P455_long.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..


sd = [2011,9,3];
ed = [2012,4,29];
figure
subplot(2,1,1)
TimeSeriesPlot('P455_timeSeries.mat',sd,ed)
subplot(2,1,2)
SnotelPlot('p455',sd,ed)
str = sprintf('P455_2012.png');
h = gcf;
cd Figures/
print(h,'-dpng',str);
cd ..

TimeSeriesYearPlot('P455_timeSeries.mat')
h = gcf;
cd Figures/
print(h,'-dpng','P455_doy.png')
cd ..
