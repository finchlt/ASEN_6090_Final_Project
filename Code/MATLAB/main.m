% Kyle Wolma
% ASEN 6090
% HW 8

%% Prep
clear all; close all;clc;
addpath('functions')

%% Load Files
fileloader();
% Feb 21
load('data\16762sp3.mat')
sp3_052=sp3;
load('data\auto0520.mat')
auto_052=auto;
load('data\p360052rinex.mat')
rinex_052=rinex;

load('data\16763sp3.mat')
sp3_053=sp3;
load('data\auto0530.mat')
auto_053=auto;
load('data\p360053rinex.mat')
rinex_053=rinex;


%% Constants
global c f0 f1 f2  usebrdc
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s
% f1 and f2
f0=10230000;    % [Hz]
f1=154*f0;      % [Hz]
f2=120*f0;      % [Hz]
% Ionosphere Free Equation
MP= @(P1,L1,L2) P1-(f1^2+f2^2)/(f1^2-f2^2)*L1+2*f2^2/(f1^2-f2^2)*L2;
 indp2=rinex_052.col.P2;
 indl1=rinex_052.col.L1;
 indl2=rinex_052.col.L2;
 
usebrdc=false;

% Sort for PRN
prn=1;
inds052=(rinex_052.data(:,3)==prn);
prnrinex_052=rinex_052.data(inds052,:);
inds053=(rinex_053.data(:,3)==prn);
prnrinex_053=rinex_053.data(inds053,:);
% Clean data 
inds052=(prnrinex_052(:,indp2)>0);
prnrinex_052=prnrinex_052(inds052,:);
inds053=(prnrinex_053(:,indp2)>0);
prnrinex_053=prnrinex_053(inds053,:);
% Calculate Multipath

 mp052=MP(prnrinex_052(:,indp2),prnrinex_052(:,indl1)*c/f1,prnrinex_052(:,indl2)*c/f2);
 mp053=MP(prnrinex_053(:,indp2),prnrinex_053(:,indl1)*c/f1,prnrinex_053(:,indl2)*c/f2);
 
 
 %% Pots
 figure
 plot(prnrinex_052(:,2)-prnrinex_052(1,2),mp052,'o')
 hold on
 plot(prnrinex_053(:,2)-prnrinex_053(1,2),mp053,'ro')
 legend('February 21, 2012','February 22, 2012')
 title('Multipath Difference')
 xlabel('Time Since 00:00 [s]')
 ylabel('MP1 [m]')
%% Path cleanup
rmpath('functions')