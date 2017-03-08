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
param.N=1000;                     % Number of noise sources
param.duration=1000;            % Source signals duration [s.]
param.temporal_sampling=0.1;    % Temporal sampling [s.]
output.F='yes';                 % Plot source power-spectrum (= FFT(auto-correlation) by Wiener-Kintchin th.)
output.signals='yes';           % Plot 5 (or less) received signals
output.setup='yes';             % Plot experimental setup
output.xcorr='yes';             % Plot cross-correlations

% Generate receivers coordinates
for i=1:param.nb_receivers
    %     param.receivers(i,:)=[0 50*(i-1) 0];
    %     param.receivers(i,:)=[0 5*(i-1) 0];
    param.receivers(i,:)=[50*(i-3) 100 0];
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

% Compute sationnary random process with Fourier method
h=param.temporal_sampling;
w=(-param.duration*(2*pi):h:param.duration*(2*pi));
n=length(w);
t=linspace(0,param.duration,n);
R=w.^2.*exp(-w.^2);
if strcmp(output.F,'yes')
    figure(2)
    semilogx(w(round((length(w))/2+1):end),R(round((length(w))/2+1):end))
    xlabel '\omega'
    ylabel 'F(\omega)'
    title 'Source Power Spectrum F(w)=w^2e^{-w^2}'
    set(gca,'FontSize',15)
end
if strcmp(output.signals,'yes')
    figure(3),hold on
end
for i=1:param.N
    W=randn(1,n);
    filter=fft(fftshift(R));
    F(i,:)=real(ifft(sqrt(filter).*fft(W)));
    if i<=5 && strcmp(output.signals,'yes')
        figure(3)
        plot(t,F(i,:))
        info{i}=sprintf('Source %d',i);
    end
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

% Compute response on each receivers
for j=1:param.nb_receivers
    for i=1:param.N
        d=norm(param.receivers(j,:)-param.sources(i,:));
        data.r{j}(i,:)=real(ifft(fft(F(i,:)).*1/(4*pi*d).*exp(1i*w*d)));
    end
    data.rtot{j}=sum(data.r{j},1);
end

% Compute empirical cross-correlation for Green's function estimation
if strcmp(output.xcorr,'yes')
    figure(4)
end
for i=1:param.nb_receivers
    [data.C(i,:),lags]=xcorr(data.rtot{1},data.rtot{i});
    lags=lags*(t(2)-t(1));
    if strcmp(output.xcorr,'yes')
        subplot(param.nb_receivers,1,i)
        plot(lags,data.C(i,:),'k')
        [~,tmp]=max(abs(data.C(i,:)));
%         xlim([lags(tmp)-10 lags(tmp)+10])
        legend(sprintf('xcorr(x_1,x_%d)',i));
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end