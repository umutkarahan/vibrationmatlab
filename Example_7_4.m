clear all 
%% Parameters
% suspended mass
m = 1;                              % mass           
% isolator
k=10000;                            % stifness  
z = 0.01; c = 2*z*sqrt(m*k);        % damping 

%% random acceleration signal
fs=2000;                            % sampling frequency
dt=1/fs;T=240;t1=0:dt:T;            % [s]
y = randn(length(t1),1);            % random acceleration signal 
N=length(y);TT=N*dt;
t=0:dt:(length(y)-1)*dt;            % time vector

%% parameters for processing the data
Na = 16;                            % number of averages
nfft=round(N/Na);                   % number of points in the DFT
noverlap=round(nfft/2);             % number of points in the overlap
df=1/(nfft*dt);                     % frequency resolution
ff = 0:df:fs/2;                     % frequency vector
w=2*pi*ff;ft=ff;

%% calculation of the acceleration of the mass
Tf=(k+j*w*c)./(k-w.^2*m+j*w*c);        % acceleration transmissibility
TTf=[Tf,fliplr(conj(Tf(1:length(Tf)-1)))];
ttf=fs*ifft(TTf);                      % IDFT of the transmissibility
x=conv(real(ttf),y)/fs;                % acceleration response   
x=x(1:length(y));                      % acceleration response 

%% calculation of the frequency domain quantities
Syy=cpsd(y,y,hann(nfft),noverlap,nfft,fs);     % PSD of shaker acceleration
Sxx=cpsd(x,x,hann(nfft),noverlap,nfft,fs);     % PSD of mass acceleration
H=tfestimate(y,x,hann(nfft),noverlap,nfft,fs); % FRF (transmissibility)
coh=mscohere(y,x,hann(nfft),noverlap,nfft,fs); % coherence

%% plot the results
%% time domain plots
figure                                          
subplot(2,1,1)
plot(t,y,'linewidth',3,'Color',[.6 .6 .6]); grid % shaker acceleration
set(gca,'fontsize',16)
axis([0,T,1.1*min(y),1.1*max(y)])
%xlabel('time (s)');
ylabel('shak. acc. (m/s^2)');
subplot(2,1,2)
plot(t,x,'linewidth',3,'Color',[.6 .6 .6]); grid  % mass acceleration
set(gca,'fontsize',16)
axis([0,T,1.1*min(x),1.1*max(x)])
xlabel('time (s)');
ylabel('mass acc. (m/s^2)');

%% Frequency domain plots
figure                                          % plot of shaker acceleration PSD
semilogx(ff,10*log10(Syy),'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,200,-50,-20])
xticks([1 10 100 1000])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('shak. acc. PSD (dB ref 1 (m^2/s^4)/Hz)');
%%
figure                                          % plot of mass acceleration PSD
semilogx(ff,10*log10(Sxx),'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,200,-80,20])
xticks([1 10 100 1000])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('mass acc. PSD (dB ref 1 (m^2/s^4)/Hz)');

figure                                             % plot of modulus of transmissibility
semilogx(ff,20*log10(abs(H)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ft,20*log10(abs(Tf)),':k','linewidth',4)
axis square; grid; axis([1,200,-60,40])
xticks([1 10 100 1000])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|Transmissibility| (dB ref unity)');

1/(2*max(abs(H)))
1/(2*max(abs(Tf)))

figure                                                    % plot of phase
semilogx(ff,180/pi*unwrap(angle(H)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ft,180/pi*unwrap(angle(Tf)),':k','linewidth',4)
xticks([1 10 100 1000])
axis square; grid; axis([1,200,-200,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                % plot of coherence
semilogx(ff,coh,'linewidth',4,'color',[0.6 0.6 0.6])
axis square; grid; axis([1,200,0,1.1])
xticks([1 10 100 1000])
yticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('coherence');