function [] = TimeSeriesYearPlot(filename)
%% Create a plot from the position
% Load data from data directory
% Inputs
% filename = 'filename.mat'
% start_date = [year,month,day]
% end_date = [year,month,day]
day = [];
year = [];
month = [];
sig =[];
soln = [];

cd ../../Data
load(filename);
cd ../Code/MATLAB

% Remove zeros
index = find(year>=2009);
day = double(day(1,index));
month = double(month(1,index));
year = double(year(1,index));
sig = sig(:,index);
soln = soln(:,index);
dates = [year',month',day'];
dec_dates = decyear(dates);

colors = [72,118,255;0,0,0;255,0,0;34,139,34;218,165,32;160,32,240;255,110,180;72,209,204;]./255;
figure
hold on
years = unique(year);
for i = 1 : length(years)
    ind{i,1} = find(year == years(i));
    soln_year{i,1} = soln(:,ind{i,1});
    sig_year{i,1} = sig(:,ind{i,1});
    dec_dates_year{i,1} = dec_dates(ind{i,1},1) - years(i);
    plot(dec_dates_year{i,1},soln_year{i,1}(3,:),'.','Color',colors(i,:),'markersize',16)
%     plot(soln_year{i,1}(3,:),dec_dates_year{i,1},'.','Color',colors(i,:))
    
    %     errorbar(dec_dates_year{i,1},soln_year{i,1}(3,:),sig_year{i,1}(3,:),'Color',colors(i,:))
end
xlabel('Day of year','fontsize',12)
ylabel('Vertical (mm)','fontsize',12)
str = sprintf('Station %s\nVerical position variation vs. day of year\nFrom: 2009 to 2012',filename(1:4));
title(str,'fontsize',12);
legend('2009','2010','2011','2012')
grid on
end