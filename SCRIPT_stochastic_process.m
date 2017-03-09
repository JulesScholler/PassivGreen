% Clear Matlab
clear
close all

% Parameters
N=1;        % Number of simulated random process
df=0.1;     % Frequency step [Hz.]
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

W=randn(n,N);
filter=repmat(R',1,N);
F=sqrt(filter).*fft(W);
for i=1:min(5,N)
    figure(2),hold on
    plot(t,real(ifft(F(:,i))))
    info{i}=sprintf('Source %d',i);
end
figure(2)
hold off
xlabel 'Time [s.]'
ylabel 'Amplitude'
title 'Noise signals'
legend(info)
clear info
set(gca,'fontsize',15)