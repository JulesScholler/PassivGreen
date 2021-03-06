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
param.N=100;                    % Number of noise sources
param.duration=10000;           % Source signals duration [s.]
param.temporal_sampling=0.1;   % Temporal sampling [s.]
output.F='no';                  % Plot source power-spectrum (= FFT(auto-correlation) by Wiener-Kintchin th.)
output.signals='no';            % Plot 5 (or less) received signals
output.setup='yes';             % Plot experimental setup
output.xcorr='yes';             % Plot cross-correlations
output.C_N='no';               % Plot C_N (cross-correlation expectations T inf)
output.C1='no';               % Plot C1 (cross-correlation expectations T inf and S inf)
tic
% Generate receivers coordinates
for i=1:param.nb_receivers
%     param.receivers(i,:)=[0 5*(i-1) 0];
    param.receivers(i,:)=[0 50*(i-1) 0];
%       param.receivers(i,:)=[50*(i-3) 100 0];
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
tau.ini=toc
tic
% Compute sationnary random process with Fourier method
h=param.temporal_sampling;
x=(-param.duration/2:h:param.duration/2); % Define grid for random process simulation
n=length(x);
t=linspace(0,param.duration,n);
f=linspace(-1/(2*h),1/(2*h),n);
w=2*pi*f;
R=x.^2.*exp(-x.^2);                 % Covariance function
if strcmp(output.F,'yes') % Plot covariance function
    figure(2)
    semilogx(x(round((length(x))/2+1):end),R(round((length(x))/2+1):end))
    xlabel '\omega'
    ylabel 'F(\omega)'
    title 'Source Power Spectrum F(w)=w^2e^{-w^2}'
    set(gca,'FontSize',15)
end
W=randn(param.N,n);         % Random white gaussian vector
filter=fft(fftshift(R));
F=sqrt(filter).*fft(W,n,2); % Generate random process with covariance R (F is in Fourier domain for now)
clear W
if strcmp(output.signals,'yes')
    figure(3),hold on
    for i=1:min(5,param.N)
        plot(t,real(ifft(F(i,:))))
        info{i}=sprintf('Source %d',i);
    end
    hold off
    xlabel 'Time [s.]'
    ylabel 'Amplitude'
    title 'Noise signals'
    legend(info)
    clear info
    set(gca,'fontsize',15)
end
tau.random_process=toc
tic
% Compute response on each receivers
Rw=(w).^2.*exp(-w.^2);
for j=1:param.nb_receivers
    C_N=zeros(1,n);
    r=zeros(1,n);
    for i=1:param.N
        d=norm(param.receivers(j,:)-param.sources(i,:)); % Distance between source and receiver
        G=1/(4*pi*d).*exp(1i*w*d);                       % Green function
        d1=norm(param.receivers(1,:)-param.sources(i,:));
        G1=1/(4*pi*d).*exp(1i*w*d1);
        % Compute response on each receivers
        r=r+real(ifft(F(i,:).*fftshift(G)));
        % Compute C_N(t,x_1,x_j) = expectation with respect to emitted signals of cross-correlation
        C_N=C_N+real(fftshift(fft(fftshift(conj(G1)).*fftshift(G).*fftshift(Rw))));
    end
    data.rtot{j}=r;
    data.C_Ntot(j,:)=C_N*(0.25*norm(param.receivers(1,:)-param.receivers(2,:)))/max(C_N);
end
clear F
tau.compute_response=toc
tic
% Compute empirical cross-correlation for Green's function estimation
for i=1:param.nb_receivers
    data.C(i,:)=real(ifftshift(ifft(fft(data.rtot{1}).*fft(fliplr(data.rtot{i})))));
%     data.C(i,:)=xcorr(data.rtot{1},data.rtot{i});
    data.C(i,:)=data.C(i,:)*(0.25*norm(param.receivers(1,:)-param.receivers(2,:)))/max(data.C(i,:));
    lags=(-n/2:(n-1)/2)*h;
%     lagx=(-n+1:n-1)*h;
    if strcmp(output.xcorr,'yes')
        figure(4),hold on
        plot(lags,data.C(i,:)+norm(param.receivers(1,:)-param.receivers(i,:)),'k')
%         plot(lagx,data.C(i,:),'k')
        legend('Empirical Cross-correlation');
    end
    if strcmp(output.C_N,'yes')
        figure(4),hold on
        plot(lags, data.C_Ntot(i,:)+norm(param.receivers(1,:)-param.receivers(i,:)), 'r')
        if strcmp(output.xcorr,'yes')
            legend('Empirical Cross-correlation','Expectation (T inf)');
        else
            legend('Expectation (T inf)');
        end
    end
    xlim([-50 300])
    set(gca,'fontsize',15)
    xlabel('Delay [s.]')
    ylabel('Distance [m.]')
end
tau.xcorr=toc