%% calculation of the FRF due to random noise
clear all 
%% Parameters
m = 1;                              % mass           
k=10000;                            % stifness  
z = 0.01; c = 2*z*sqrt(m*k);        % damping 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % natural frequency
fn=wn/(2*pi);                       % natural frequency
Tn=1/fn;                            % natural period 
SNRf=1000;SNRx=1000;                % added noise SNRs

%% Random force signal
fs=2000;                           % Sampling frequency
dt=1/fs;T=120;t=0:dt:T;            % time vector
f=randn(length(t),1);              % random signal 
N=length(f);

%% Calculation of displacement response
% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % IRF
% Convolution
xc=conv(h,f)/fs;                    % displacement response   
xc=xc(1:length(f));

%% Frequency domain calculations
fwn=awgn(f,SNRf,'measured','dB');                  % add random noise
xwn=awgn(xc,SNRx,'measured','dB');                 % add random noise

Na = 32;                            % number of averages
nfft=round(N/Na);                  % number of points in the DFT
noverlap=round(nfft/2);            % number of points in the overlap
Sff=cpsd(fwn,fwn,hann(nfft),noverlap,nfft,fs); % PSD
Sxx=cpsd(xwn,xwn,hann(nfft),noverlap,nfft,fs); % PSD
Sfx=cpsd(xwn,fwn,hann(nfft),noverlap,nfft,fs); % PSD
Tfx=tfestimate(fwn,xwn,hann(nfft),noverlap,nfft,fs); % FRF
Coh=mscohere(fwn,xwn,hann(nfft),noverlap,nfft,fs); % coherence
df=1/(nfft*dt);                     % frequency resolution
ff = 0:df:fs/2;                     % frequency vector
       
%% theory
df=0.001;ft=0:df:fs-df;             % frequency vector
w=2*pi*ft;
Ht=1./(k-w.^2*m+j*w*c);             % theoretical FRF

%% Time domain plots
figure
subplot(2,1,1)
plot(t,fwn,'linewidth',3,'Color',[.6 .6 .6]),grid
set(gca,'fontsize',16)
axis([0,T,1.1*min(fwn),1.1*max(fwn)])
xlabel('time (s)');
ylabel('force (N)');
subplot(2,1,2)
plot(t,xwn,'linewidth',3,'Color',[.6 .6 .6]),grid
set(gca,'fontsize',16)
axis([0,T,1.1*min(xwn),1.1*max(xwn)])
xlabel('time (s)');
ylabel('displacement (m)');

%% Frequency domain plots
figure                                          % plot of force PSD
semilogx(ff,10*log10(Sff),'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,1000,-60,-20])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('force PSD (dB ref 1 N^2/Hz)');

figure                                          % plot of displacement PSD
semilogx(ff,10*log10(Sxx),'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,1000,-180,-60])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('displ. PSD (dB ref 1 m^2/Hz)');

figure                                               % plot of modulus of FRF
semilogx(ff,20*log10(abs(Tfx)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ft,20*log10(abs(Ht)),':k','linewidth',4)
axis square; grid; axis([1,1000,-140,-40])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('displ./force (dB ref 1 m/N)');

figure                                              % plot of phase
semilogx(ff,180/pi*angle(Tfx),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ft,180/pi*unwrap(angle(Ht)),':k','linewidth',4)
axis square; grid; axis([1,1000,-200,0])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                % plot of coherence              
semilogx(ff,Coh,'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,1000,0,1])
xticks([1 10 100 1000])
yticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('coherence');
