close all
clear all

eps_i=10;
eps_a=1;
eps_w=0.5;

u = 2; % form number

P_i=0.5;
P_a=0.4;
P_w=0.1;

tot_vol=1; % m^3
mass_i=1; % kg

i=1;

for p_i = 0.2:0.01:0.8
    p_a = 0.9-p_i;
    
    m = P_i*((eps_i-1)/(eps_i+u)) + P_a*((eps_a-1)/(eps_a+u)) + ...
        P_w*((eps_w-1)/(eps_w+u));

    eps_s(i) = (m*u-1)/(m-1);
    i=i+1;
end

plot(0.2:0.01:0.8,eps_s)
