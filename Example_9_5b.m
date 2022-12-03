%clear all

%% parameters
% structure
m=1;                                         % mass
k=1e4;                                       % stiffness
wn=sqrt(k/m);                                % natural frequency
% absorber
ma=0.05;                                     % mass of absorber
mu=ma/m;                                     % mass ratio
wa=wn/(1+mu);                                % absorber natural frequency
ka=wa^2*ma;                                  % stiffness of the absorber
za=sqrt(3/8*mu/(1+mu)^3);                    % damping ratio
ca=2*za*sqrt(ma*ka);                         % damping coefficient
                                   
%shaker
ms = 0.1;                                     % mass  
ws=2*pi*10;                                   % natural frequency (rad/s)
ks=ws^2*ms;                                   % stiffness
zs=0.1; cs=2*zs*sqrt(ms*ks);                  % damping

%% direct FRF calculation
fs=1000;df=0.01;dt=1/fs;                     % frequency parameters
f=0:df:fs/2;w=2*pi*f;                        % frequency vector
Ka=ma*(ka+j*w*ca)./(ka-w.^2*ma+j*w*ca);      % apparent mass
                                             
%% IRF
Hd=[Ka fliplr(conj(Ka))];                    % form the double-sided spectrum  
Ht=Hd(1:length(Hd)-1);                       % set the length of the FRF
h=fs*real(ifft(Ht));                         % calculation of the IRF

%% input and output
T=250;t=0:dt:T;                              % signal duration; time vector                  
a = randn(1,length(t));                      % acceleration random signal 
a=a-mean(a);                                 % set the mean to zero
N=length(a);
fe = conv(h,a)/fs;                            % force signal   
fe = fe(1:length(a));                         % force signal   

%% frequency domain calculations
Na = 8;                                           % number of averages
nfft=round(N/Na);                                 % number of points in the DFT
noverlap=round(nfft/2);                           % number of points in the overlap
Sff=cpsd(fe,fe,hann(nfft),noverlap,nfft,fs);      % PSD of force signal
Saa=cpsd(a,a,hann(nfft),noverlap,nfft,fs);        % PSD of acceleration signal
Kae=tfestimate(a,fe,hann(nfft),noverlap,nfft,fs); % apparent mass
CoHe=mscohere(a,fe,hann(nfft),noverlap,nfft,fs);  % coherence
dff=1/(nfft*dt);                                  % frequency resolution
ff = 0:dff:fs/2;                                  % frequency vector
%% plot the results
figure                          % force time history
plot(t,fe,'linewidth',2,'color',[0.6 0.6 0.6]) 
axis square,axis([0,250,-0.2,0.2]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('force (N)');

figure                         % acceleration time history
plot(t,a,'linewidth',2,'color',[0.6 0.6 0.6])  
axis square,axis([0,250,-6,6]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('acceleration (m/s^2)');

figure                                               % force PSD
semilogx(ff,10*log10(abs(Sff)),'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,-100,-30]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('force PSD (dB ref 1 N^2/Hz)');

figure                                               % acceleration PSD
semilogx(ff,10*log10(abs(Saa)),'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,-50,-10]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('acc. PSD (dB ref 1 m^2/s^4Hz)');

figure                                               % apparent mass
semilogx(ff,20*log10(abs(Kae)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,20*log10(abs(Ka)),'--k','linewidth',4)
axis square,axis([1,1000,-70,-10]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|apparent mass| (dB ref 1 kg)');

figure                                                      % phase
semilogx(ff,180/pi*unwrap(angle(Kae)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,180/pi*unwrap(angle(Ka)),'--k','linewidth',4)
axis square,axis([1,1000,-200,0]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                 % coherence
semilogx(ff,CoHe,'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,0,1]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('coherence');

figure                                         % Nyquist plot
plot(real(Kae),imag(Kae),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(real(Ka),imag(Ka),'--k','linewidth',4)
axis square,axis([-0.15,0.15,-0.3,0]),grid
set(gca,'fontsize',16)
xlabel('real \{apparent mass\} (kg)');
ylabel('imag \{apparent mass\} (kg)');
