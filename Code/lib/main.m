%%  main.m
clear all; clc; close all;
%addpath ..\code4students; 
c = 2.99792458e8;
%% 1) Read header file, extract the receiver position
[~,userPos,~] = read_rinex_header('namb286x.07o');
userPos = [userPos; 0];
%% 2) Read the first set of observation records(One epoch, all satellites).
% For each satellite, calculate the following values
%   a. the ionospheric free pseudorange observable (using P1 and P2)
%   b. the expected range
%   c. the elevation and azimuth
%   d. the satellite clock correction (using a0 and a1)
%   e. the relativistic correction

% Extract the ephemeris, using yesterday's
[eph,~] = read_GPSbroadcast('auto2860.07n');
% Extract 100 rows of observations
obs = read_rinex_obs('namb286x.07o',1:32,100);
[lat, lon, alt] = wgsxyz2lla(userPos');
% Cut the first epoch
index = find(obs.data(:,2) == obs.data(1,2));
satNum = length(index);
obsVar = obs.data(index,:);

satPos = zeros(satNum,3);
satVel = zeros(satNum,3);
Pseudorange_ionofree = zeros(1,satNum);
Georange_expected = zeros(1,satNum);
azi = zeros(1,satNum);
ele = zeros(1,satNum);
dtr = zeros(1,satNum);
dtsv = zeros(1,satNum);
dx = 100*ones(4,1);
t_mear = 0;
while  norm(dx(1:4)) > 1e-4
    for i = 1:satNum
        %   Calculate the ionospheric free pseudorange observalbe
        f1 = 1575.42e6;
        f2 = 1227.6e6;
        Pseudorange_ionofree(i) = f1^2/(f1^2 - f2^2)*obsVar(i,6) - f2^2/(f1^2 - f2^2)*obsVar(i,7);
        %   Calculate the expected range
        prn = obsVar(i,3);
        t = obsVar(i,1:2);
        t(2) = t(2) + userPos(4)/c;
        [Georange_expected(i),tsv,~] = compute_range(eph, prn, t,userPos(1:3)'); 
        %   Calculate the elevation and azimuth, at the transmission time
        [h, satPos(i,:), satVel(i,:),currEph] = broadcast2xv1(eph,[t(1) tsv],prn);
        [azi(i), ele(i),~] = ecef2azelrange2(satPos(i,:),userPos,lat,lon);
        %   Calculate the relativistic correction
        dtr(i) = - 2*satPos(i,:)*satVel(i,:)'/(c^2);
        %   Calculate the satellite clokc error
        toc = currEph(17);
        af0 = currEph(21);
        af1 = currEph(22);
        af2 = currEph(23);
        dt = tsv - toc;
        if dt > 302400
            dt = dt - 604800;
        end
        if dt < -302400
            dt = dt + 604800;
        end
        dtsv(i) = (af0 + af1*dt + af2*dt^2 + dtr(i))*c;    
    end

    %% Construct the A matrix
    A = zeros(satNum, 4);
    x0 = userPos(1);
    y0 = userPos(2);
    z0 = userPos(3);
    for i = 1: satNum
         A(i,1) = -(satPos(i,1) - x0)/Georange_expected(i);
         A(i,2) = -(satPos(i,2) - y0)/Georange_expected(i);
         A(i,3) = -(satPos(i,3) - z0)/Georange_expected(i);
         A(i,4) = 1;
    end
    %% Compute dy
    dy = Pseudorange_ionofree' - Georange_expected' + dtsv'- userPos(4);
    %% Least Square
    dx = A\dy;
    userPos_cal = userPos + dx(1:4);
    t_mear = userPos_cal(4)/c;
    e = dy - A*dx;
    userPos = userPos_cal;
    [wlat, wlon, walt] = wgsxyz2lla(userPos_cal);
    %enu = wgslla2enu(wlat, wlon, walt, lat, lon, alt);
end

figure(1);
plot(ele,dy,'.r');

