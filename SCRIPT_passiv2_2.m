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
param.N=10;                     % Number of noise sources
param.duration=1000;            % Source signals duration [s.]
param.temporal_sampling=0.05;    % Temporal sampling [s.]
output.F='no';                 % Plot source power-spectrum (= FFT(auto-correlation) by Wiener-Kintchin th.)
output.signals='no';           % Plot 5 (or less) received signals
output.setup='no';             % Plot experimental setup
output.xcorr='yes';             % Plot cross-correlations
output.C_N='yes';               % Plot C_N (cross-correlation expectations)
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
F=sqrt(filter).*fft(W,n,2); % F^(i,w) = realization of the random process n^(x_i,w)
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
% Compute response on each receiver
for j=1:param.nb_receivers
    for s=1:param.N
        d=norm(param.receivers(j,:)-param.sources(s,:))/param.duration*n;
        data.r{j}(s,:)=real(ifft(F(s,:).*1/(4*pi*d).*exp(1i*f*pi*h*d)));
        data.G{j}(s,:)=1/(4*pi*d).*exp(1i*f*pi*h*d);
    end
    data.rtot{j}=sum(data.r{j},1);
end
tau.compute_response=toc;
tic
% Compute empirical cross-correlation for Green's function estimation
for j=1:param.nb_receivers
    data.C(j,:)=real(ifftshift(ifft(fft(data.rtot{1}).*fft(fliplr(data.rtot{j})))));
    lags=(-length(data.C(j,:))/2:(length(data.C(j,:))-1)/2)*h;
    if strcmp(output.xcorr,'yes')
        figure(4)
        subplot(param.nb_receivers,1,j)
        plot(lags, data.C(j,:),'k')
        [~,tmp]=max(abs(data.C(j,:)));
        %         xlim([lags(tmp)-100 lags(tmp)+100])
        legend(sprintf('xcorr(x_1,x_%d)',j));
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end
tau.xcorr=toc;
tic
% Compute C_N(t,x_1,x_j) = expectation with respect to emitted signals of cross-correlation
for j=1:param.nb_receivers
    for s = 1 : param.N
        data.C_N{j}(s,:) = real(fftshift(fft(conj(data.G{1}(s,:)).*data.G{j}(s,:).*R)));
    end
    data.C_Ntot(j,:) = sum(data.C_N{j}, 1);
    lags=(-n/2:(n-1)/2)*h;
    if strcmp(output.C_N,'yes')
        figure(5)
        subplot(param.nb_receivers,1,j)
        plot(lags, data.C_Ntot(j,:), 'k')
        legend(sprintf('C_N(t,x_1,x_%d)',j));
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end
tau.C_N=toc;

tic
% Compute C1(t,x_1,x_j) expectation with respect to the distribution of emitted signals 
% and the source positions of cross-correlation
d1 = @(y1, y2, y3) norm(param.receivers(1,:)-[y1, y2, y3])/param.duration*n;
%G1 = @(y1, y2, y3) 1/(4*pi*d1(y1, y2, y3)).*exp(1i*f*pi*h*d1(y1, y2, y3)); %G(w, x_1, y) (n-vector)
%K = @(y1, y2, y3) mvnpdf([y1, y2, y3], param.mu, param.sigma);
K = @(y1, y2, y3) normpdf(y1, param.mu(1), param.sigma(1)) .* normpdf(y2, ...
    param.mu(2), param.sigma(2)) .* normpdf(y3, param.mu(3), param.sigma(3));

for j=1:param.nb_receivers
    d = @(y1, y2, y3) norm(param.receivers(j,:)-[y1, y2, y3])/param.duration*n;
    temp = zeros(1,n);
    for i = 1:n
        G1 = @(y1, y2, y3) 1/(4*pi*d1(y1, y2, y3)).*exp(1i*f(i)*pi*h*d1(y1, y2, y3)); %G(w_i, x_1, y)
        G = @(y1, y2, y3) 1/(4*pi*d(y1, y2, y3)).* exp(1i*f(i)*pi*h*d(y1, y2, y3)); %G(w_i, x_j, y)
        func = @(y1, y2, y3) K(y1, y2, y3) * conj(G1(y1, y2, y3)) .* G(y1, y2, y3);
        temp(i) = integral3(func,-Inf,Inf,-Inf,Inf,-Inf,Inf);
    end
    data.C1(j,:) = real(fftshift(fft(temp.*R)));p

    lags=(-n/2:(n-1)/2)*h;
    if strcmp(output.C1,'yes')
        figure(6)
        subplot(param.nb_receivers,1,j)
        plot(lags, data.C1(j,:), 'k')
        legend(sprintf('C1(t,x_1,x_%d)',j));
        set(gca,'fontsize',15)
        xlabel('Delay [s.]')
        ylabel('Ampl.')
    end
end
tau.C1=toc