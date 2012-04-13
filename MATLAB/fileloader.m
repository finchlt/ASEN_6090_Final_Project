% Kyle Wolma
% ASEN 6090
% Final Project
% File loader

function fileloader()
% Sp3
if(exist('Data\16762sp3.mat'))
else
    sp3=read_sp3('Data\igs16762.sp3');
    save Data\16762sp3.mat
    clear all;
end
if(exist('Data\16763sp3.mat'))
else
    sp3=read_sp3('Data\igs16763.sp3');
    save Data\16763sp3.mat
    clear all;
end

% Navigation Files
if(exist('Data\auto0520.mat'))
else
    auto=read_GPSbroadcast('Data\auto0520.12n');
    save Data\auto0520.mat
    clear all;
end
if(exist('Data\auto0530.mat'))
else
    auto=read_GPSbroadcast('Data\auto0530.12n');
    save Data\auto0530.mat
    clear all;
end
% Rinex Observations
if(exist('Data\p360052rinex.mat'))
else
    rinex=read_rinex_obs('Data\p3600520.12o');
    save Data\p360052rinex.mat
    clear all;
end

if(exist('Data\p360053rinex.mat'))
else
    rinex=read_rinex_obs('Data\p3600530.12o');
    save Data\p360053rinex.mat
    clear all;
end
