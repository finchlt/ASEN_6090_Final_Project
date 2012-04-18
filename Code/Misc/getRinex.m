%%% Get rinex data 
%%% Praveen Vikram

clear all
close all
% clc;

cur_folder = [pwd '/'];
lib_folder = [pwd '/../lib'];
data_folder = [pwd '/../../Data'];
py_folder = [pwd '/../py'];

addpath(lib_folder);

station='p360';
year='2012';
doy='052';

cmd=[py_folder '/getFiles.py --station=' station ' --year=' year ' --doy=' doy];
disp(cmd)

system('pwd')
system(cmd);
