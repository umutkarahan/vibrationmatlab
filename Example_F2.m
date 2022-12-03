
clear all 
%% parameters
m = 10;                              % mass           
k=1000;                              % stiffness  
z = 0.5; c = 2*z*sqrt(m*k);          % damping 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;     % undamped and damped natural frequencies

%% time and frequency parameters
T=10;                                % time window
fs=10;                               % sampling frequency
dt=1/fs; t=0:dt:T;                   % time resolution and time vector
df=1/T; f=0:df:fs;                   % frequency resolution and frequency vector
N=length(t);                             

%% impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);  % IRF

%% calculation of DFT
H=dt*fft(h);                         % calculation of receptance FRF

%% calculation of aliased response
for p=1:1000                         % calculation of the frequency range that is included in the sum of aliases for the theoretical FRF
    f1=(p-1)*fs:df:p*fs;
    w1=2*pi*f1;
    HP(p,:)=1./(k-w1.^2*m+j*w1*c);   % aliased FRF for +ve frequencies
    f2=-p*fs:df:-(p-1)*fs;
    w2=2*pi*f2;
    HM(p,:)=1./(k-w2.^2*m+j*w2*c);   % aliased FRF for -ve frequencies
end
MP=sum(HP); MS=(sum(HM));            % summing all aliases
HT1=MP+MS;                           % total aliased FRF

HT=HT1(1:(N+1)/2);
HA=[HT fliplr(conj(HT))];            % double sided spectrum
HA1=HA(1:length(f));

%% plot of results
figure                                    % modulus
plot(f,abs(HA1),'k','linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,abs(H),'ok','linewidth',2,'markersize',8)
xticks([0 2 4 6 8 10])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (m/N)');
grid;axis square

figure                                             % phase
plot(f,180/pi*angle(HA1),'k','linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,180/pi*angle(H),'ok','linewidth',2,'markersize',8)
xticks([0 2 4 6 8 10])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
grid;axis square