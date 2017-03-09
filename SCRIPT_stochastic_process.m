% Clear Matlab
clear
close all

% Parameters
N=1;        % Number of simulated random process
df=0.01;     % Frequency step [Hz.]
fmax=10;    % Maximum frequency [Hz.]

w=(-fmax:df:fmax);
n=length(w);
t=0:1/df:(length(w)-1)/df;
R=w.^2.*exp(-w.^2);

% Plot F
figure(1)
semilogx(w(round((length(w))/2+1):end),R(round((length(w))/2+1):end))
xlabel '\omega'
ylabel 'F(\omega)'
title 'Source Power Spectrum F(w)=w^2e^{-w^2}'
set(gca,'FontSize',15)

W=randn(n,1);
filter=fft(fftshift(R))';
F=sqrt(filter).*fft(W);

% Plot random process
figure(2)
plot(t,real(ifft(F)))
xlabel 'Time [s.]'
ylabel 'Amplitude'
title 'Noise signals'
set(gca,'fontsize',15)