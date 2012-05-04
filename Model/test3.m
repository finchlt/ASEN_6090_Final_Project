close all
clear all

f = 1.57542e9;
c = 299792500;
lambda = c/f;

% T = -20.55;
T = 0;
rho_d = 0.24;
eps_re = 1 + 2*rho_d;
eps_im = eps_re * 1.59 * 1e6 * (0.52*rho_d + 0.62*rho_d^2)/(7 + 1.7*rho_d + 0.7*rho_d^2)...
	*(f^-1 + 1.23*1e-14*sqrt(f))*exp(0.036*T);

eps_j = eps_re - j*eps_im;

%%
%eps_j=1.48-j*2.71e-4; % complex dielectric coeff of snow

%eps_j=1.48-j*2.71e-4;
%eps_j=2.93-j*8e-4; % for ice at GPS freq

% d is the snow depth on antena (m)

omega=2*pi*f; % angular freq
k = omega*sqrt(eps_j)/c;
v = c/sqrt(eps_j); % wave velocity in snow


i=1;
d_strt = 0.0;
d_step = 0.001;
d_end = 2;

for d = d_strt:d_step:d_end
    t(i) = d/abs(v);
    E(i) = exp(j*( real(k)*d - omega*t(i) ))*exp(imag(k)*d);
    A(i) = 1 - abs(E(i));
    i=i+1;
end

%plot(abs(E),d_end:-d_step:d_strt)
fontSize=14;

plot(d_strt:d_step:d_end,abs(E))
xlim([d_strt d_end])
title('Modeled EM wave Attenuation','FontSize',fontSize)
ylabel('Attenuation','FontSize',fontSize)
xlabel('Snow thickness','FontSize',fontSize)
set(gca,'FontSize',fontSize)
grid on
text(1.2,0.9998,'\rho_d=0.24','FontSize',fontSize)
text(1.2,0.9997,'T=0C','FontSize',fontSize)
text(1.2,0.9996,'\epsilon=1.48-j5.69x10^{-5}','FontSize',fontSize)

