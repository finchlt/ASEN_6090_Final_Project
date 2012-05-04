close all
clear all

f = 1.57542e9;
c = 299792500;
lambda = c/f;

theta=pi/2; % angle of incidence (assume directly overhead)
d=0.5; % snow depth on antena (m)
eps_j=1.48-j*2.71e-4; % complex dielectric coeff of snow

i=1;
for d = 0.1:0.01:2
    
    phi_j = (2*pi/lambda)*d*sqrt(eps_j- (cos(theta))^2 );

    Z_hj = sin(theta)/sqrt(eps_j- (cos(theta))^2);
    Z_vj = 1/Z_hj;
    Z_h3 = 1;
    Z_v3 = 1;

    W_h = (Z_h3 + j*Z_hj*tan(phi_j))/(1+j*(Z_h3/Z_hj)*tan(phi_j)); 
    W_v = (Z_v3 + j*Z_vj*tan(phi_j))/(1+j*(Z_v3/Z_vj)*tan(phi_j)); 

    r_h = (W_h-1)/(W_h+1);
    r_v = (W_v-1)/(W_v+1);

    P(i) = abs(1 + (r_h-r_v)/2)^2;
    i=i+1;
end

plot(0.1:0.01:2,P,'o-')