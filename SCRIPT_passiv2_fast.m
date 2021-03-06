% Script for passiv2
%
%
% Jules Scholler - Mar. 2017

% Reset Matlab
close all
clear

% Parameters
param.nb_receivers=5;           % Number of receivers
param.sigma=[100 50 100];       % Sources std position
param.mu=[0 -200 0];            % Sources mean position
param.N=100;                     % Number of noise sources
param.duration=10000;            % Source signals duration [s.]
param.temporal_sampling=0.05;    % Temporal sampling [s.]
output.F='no';                 % Plot source power-spectrum (= FFT(auto-correlation) by Wiener-Kintchin th.)
output.signals='no';           % Plot 5 (or less) received signals
output.setup='no';             % Plot experimental setup
output.xcorr='no';             % Plot cross-correlations
output.C_N='no';               % Plot C_N (cross-correlation expectations)
output.C1='no';               % Plot C1 (cross-correlation expectations)
tic
% Generate receivers coordinates
for i=1:param.nb_receivers
            param.receivers(i,:)=[0 5*(i-1) 0];
%     param.receivers(i,:)=[0 50*(i-1) 0];
%         param.receivers(i,:)=[50*(i-3) 100 0];
    C(i,:)=[0 0 1]; % Receivers are blue
end
if strcmp(output.setup,'yes')
    figure(1), hold on, grid on
    scatter3(param.receivers(:,1),param.receivers(:,2),param.receivers(:,3),10,C);
    xlabel 'x'
    ylabel 'y'
    zlabel 'z'
    title 'Experimental setup'
    clear C
end
% Compute and plot sources position
for i=1:param.N
    for j=1:3
        param.sources(i,j)=param.sigma(j)*randn(1)+param.mu(j);
    end
    C(i,:)=[1 0 0]; % Sources are red
end
if strcmp(output.setup,'yes')
    figure(1)
    scatter3(param.sources(:,1),param.sources(:,2),param.sources(:,3),5,C);
    legend('Receivers','Sources')
    clear C
    set(gca,'FontSize',15)
end
tau.ini=toc;
tic
% Compute sationnary random process with Fourier method
h=param.temporal_sampling;
w=(-param.duration:h:param.duration);
n=length(w);
t=(w-min(w))/2;
f=linspace(-1/(2*h),1/(2*h),n);
R=w.^2.*exp(-w.^2);
if strcmp(output.F,'yes')
    figure(2)
    semilogx(w(round((length(w))/2+1):end),R(round((length(w))/2+1):end))
    xlabel '\omega'
    ylabel 'F(\omega)'
    title 'Source Power Spectrum F(w)=w^2e^{-w^2}'
    set(gca,'FontSize',15)
end
W=randn(param.N,n);
filter=fft(fftshift(R));
F=sqrt(filter).*fft(W,n,2);
if strcmp(output.signals,'yes')
    figure(3),hold on
    for i=1:min(5,param.N)
        plot(t,real(ifft(F(i,:))))
        info{i}=sprintf('Source %d',i);
    end
    hold off
end
if strcmp(output.signals,'yes')
    figure(3)
    hold off
    xlabel 'Time [s.]'
    ylabel 'Amplitude'
    title 'Noise signals'
    legend(info)
    clear info
    set(gca,'fontsize',15)
end
tau.random_process=toc;
tic
% Compute response on each receivers
r=zeros(n,param.N,param.nb_receivers);
GreenF=r;
rtot=zeros(n,param.nb_receivers);
for j=1:param.nb_receivers
    for i=1:param.N
        d=norm(param.receivers(j,:)-param.sources(i,:)); % Distance between source and receiver
        r(:,i,j)=real(ifft(F(i,:).*1/(4*pi*d).*fftshift(exp(1i*w*d))));  % Response after propagation
        GreenF(:,i,j)=1/(4*pi*d).*exp(1i*w*d);
    end
    rtot(:,j)=sum(r(:,:,j),2);
end
tau.compute_response=toc;
tic
% Compute empirical cross-correlation for Green's function estimation
if strcmp(output.xcorr,'yes')
    figure(4)
end
C=zeros(n,param.nb_receivers);
for i=1:param.nb_receivers
    C(:,i)=real(ifftshift(ifft(fft(rtot(:,1)).*fft(fliplr(rtot(:,i))))));
    C(:,i)=C(:,i)/max(C(:,i));
    lags=(-n/2:(n-1)/2)*h;
    if strcmp(output.xcorr,'yes')
        subplot(param.nb_receivers,1,i),hold on
        plot(lags,C(:,i),'k')
        [~,tmp]=max(abs(C(:,i)));
%         xlim([lags(tmp)-100 lags(tmp)+100])
        legend(sprintf('xcorr(x_1,x_%d)',i));
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end
tau.xcorr=toc;

% Compute C_N(t,x_1,x_j) = expectation with respect to emitted signals of cross-correlation
tic
Rw=((w).^2.*exp(-w.^2))';
C_N=zeros(n,param.N,param.nb_receivers);
C_Ntot=zeros(n,param.nb_receivers);
for j=1:param.nb_receivers
    for s=1:param.N
        C_N(:,s,j) = real(fftshift(fft(fftshift(conj(GreenF(:,s,1))).*fftshift(GreenF(:,s,j)).*fftshift(Rw))));
    end
    C_Ntot(:,j) = sum(C_N(:,:,j),2);
    C_Ntot(:,j)=C_Ntot(:,j)/max(C_Ntot(:,j));
    if strcmp(output.C_N,'yes')
        figure(4)
        subplot(param.nb_receivers,1,j)
        plot(lags, C_Ntot(:,j), 'r')
        if strcmp(output.xcorr,'yes')
            legend(sprintf('Xcorr(t,x_1,x_%d)',j),sprintf('C_N(t,x_1,x_%d)',j));
        else
        legend(sprintf('C_N(t,x_1,x_%d)',j));
        end
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end
tau.C_N=toc