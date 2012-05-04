%% Load Data

clc;
clear all; 
close all;

cd ..
cd ..
cd Data

feb21 = 'p3600520.12o';
feb22 = 'p3600530.12o';
sp3_21 = 'igs16762.sp3';
sp3_22 = 'igs16763.sp3';

[ ~, rec_pos_p360, ~] = read_rinex_header(feb21);
[  p360_feb21 ] = read_rinex_obs(feb21);

[ ~, rec_pos_p360, ~] = read_rinex_header(feb21);
[  p360_feb22 ] = read_rinex_obs(feb22);

 sp3_21 = read_sp3(sp3_21);
 sp3_22 = read_sp3(sp3_22);

clear feb*

cd ..
cd Code
cd MATLAB
%% Save
save snr_data.mat