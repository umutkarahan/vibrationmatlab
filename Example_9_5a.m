clear all

%% parameters
%structure
m=1;                                         % mass
k=1e4;                                       % stiffness
wn=sqrt(k/m);                                % natural frequency
z=0.01;                                      % damping ratio
c=2*z*sqrt(k*m);                             % damping coefficient
                                    
%shaker
ms = 0.1;                                     % mass  
ws=2*pi*10;                                   % natural frequency (rad/s)
ks=ws^2*ms;                                   % stiffness
zs=0.1; cs=2*zs*sqrt(ms*ks);                  % damping

%% direct FRF calculation
fs=1000;df=0.01;dt=1/fs;                      % frequency parameters
f=0:df:fs/2; w=2*pi*f;                        % frequency vector
H=1./(k-w.^2*m+j*w*c);                        % receptance
Ha=-w.^2.*H;                                  % accelerance

%% IRF
Hd=[Ha fliplr(conj(Ha))];                      % form the double-sided spectrum  
Ht=Hd(1:length(Hd)-1);                       % set the length of the FRF
h=fs*ifft(Ht);                               % calculation of the IRF

%% input and output
T=250;t=0:dt:T;                              % signal duration; time vector                  
fe = randn(1,length(t));                     % random signal 
fe=fe-mean(fe);                              % set the mean to zero
K=1./H;                                      % dynamic stiffness of structure
Ks=ks-w.^2*ms+j*w*cs;                        % dynamic stiffness of shaker
Fe=K./(K+Ks);                                % amount of force applied to structure
G=[Fe,fliplr(conj(Fe(1:length(Fe)-1)))];     % double sided spectrum
g=fs*ifft(G);                                % IRF to determine force applied to structure
fes = conv(real(g),fe)/fs;                   % force applied to structure   
fes = fes(1:length(fe));                     % force applied to structure 
N=length(fes);
a = conv(h,fes)/fs;                          % acceleration response   
a = a(1:length(fe));                         % acceleration response   

%% frequency domain calculations
Na = 8;                                            % number of averages
nfft=round(N/Na);                                  % number of points in the DFT
noverlap=round(nfft/2);                            % number of points in the overlap
Sff=cpsd(fes,fes,hann(nfft),noverlap,nfft,fs);     % PSD of force applied to the structure
Sxx=cpsd(a,a,hann(nfft),noverlap,nfft,fs);         % PSD of acceleration response
He=tfestimate(fes,a,hann(nfft),noverlap,nfft,fs);  % FRF
CoHe=mscohere(fes,a,hann(nfft),noverlap,nfft,fs);  % coherence
dff=1/(nfft*dt);                                   % frequency resolution
ff = 0:dff:fs/2;                                   % frequency vector
%% plot the results
figure                           % force time history
plot(t,fes,'linewidth',2,'color',[0.6 0.6 0.6]) 
axis square,axis([0,250,-6,6]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('force (N)');

figure                         % acceleration time history
plot(t,a,'linewidth',2,'color',[0.6 0.6 0.6])  
axis square,axis([0,250,-8,8]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('acceleration (m/s^2)');

figure                                               % force PSD
semilogx(ff,10*log10(abs(Sff)),'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,-50,-10]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('force PSD (dB ref 1 N^2/Hz)');

figure                                               % acceleration PSD
semilogx(ff,10*log10(abs(Sxx)),'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,-80,20]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('acc. PSD (dB ref 1 m^2/s^4Hz)');

figure                                              % FRF modulus
semilogx(ff,20*log10(abs(He)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(f,20*log10(abs(Ha)),'--k','linewidth',4)
axis square,axis([1,1000,-40,40]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|accelerance| (dB ref 1 m/Ns^2)');

figure                                                     % phase
semilogx(ff,180/pi*unwrap(angle(He)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,180/pi*unwrap(angle(Ha)),'--k','linewidth',4)
axis square,axis([1,1000,0,200]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                 % coherence
semilogx(ff,CoHe,'linewidth',4,'color',[0.5 0.5 0.5])
axis square,axis([1,1000,0,1]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('coherence');

figure                                       % Nyquist plot
plot(real(He),imag(He),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(real(Ha),imag(Ha),'--k','linewidth',4)
axis square,axis([-30,30,0,60]),grid
set(gca,'fontsize',16)
xlabel('real \{accelerance\} (m/Ns^2)');
ylabel('imag \{accelerance\} (m/Ns^2)');
