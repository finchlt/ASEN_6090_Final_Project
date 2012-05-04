function [] = TimeSeriesPlot(filename,start_date,end_date)
%% Create a plot from the position
global fsize asize
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
index = find(year>0);
day = double(day(1,index));
month = double(month(1,index));
year = double(year(1,index));
sig = sig(:,index);
soln = soln(:,index);

years = unique(year);
for i = 1 : length(years)
    ind{i,1} = find(year == years(i));
end
   

dates = [year',month',day'];
dec_dates = decyear(dates);
start_date_dec = decyear(start_date);
end_date_dec = decyear(end_date);

% If no input for start_date or end_date 
if nargin < 3
    end_date = dates(end,:);
    end_date_dec = dec_dates(end,1);
end
if nargin < 2
    start_date = dates(1,:);
    start_date_dec = dec_dates(1,1);
end

% Parse data to only include between start and end date
index = find(dec_dates >= start_date_dec & dec_dates <= end_date_dec);

% SerialDates = SerialDates(index,1);
% dates = dates(index);
dates_dec = dec_dates(index);
dates = dates(index,:);
sig = sig(:,index);
soln = soln(:,index);

SerialDates = datenum(dates(:,1),dates(:,2),dates(:,3));
SerialStart = datenum(start_date(1),start_date(2),start_date(3));
SerialEnd = datenum(end_date(1),end_date(2),end_date(3));

% Plot Data
% errorbar(dates_dec,soln(3,:),sig(3,:),'*');
plot(dates_dec,soln(3,:),'.','markersize',16);
% datetick('x','yyyy.dd');
set(gca,'FontSize',asize)
xlabel('Year','fontsize',fsize)
ylabel('Vertical (mm)','fontsize',fsize)
datestr_start = datestr(SerialDates(1),'yyyy-mm-dd');
datestr_end = datestr(SerialDates(end),'yyyy-mm-dd');
str = sprintf('Station %s\nPosition time series\nFrom: %s to %s',filename(1:4),datestr_start,datestr_end);
title(str,'fontsize',fsize);
grid on
axis fill
h = gcf;
cd Figures/
% datestr_start = datestr(SerialDates(1),'yyyy-mm-dd');
% datestr_end = datestr(SerialDates(end),'yyyy-mm-dd');
cd ..

end

