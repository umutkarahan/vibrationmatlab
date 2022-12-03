
clear all

%% time and frequency parameters
fs=200;                            % sampling frequency 
T = 5;dt=1/fs;t=0:dt:T;            % signal duration; time resolution; time vector                  
N=length(t);A=1;                   % number of points; amplitude of excitation

%% random excitation
x = A*randn(length(t),1);          % random signal 
x=x-mean(x);                       % set the mean to zero

%% calculation of PSD
Na = 4;                            % number of averages
nfft=round(N/Na);                  % number of points in the DFT
noverlap=round(nfft/2);            % number of points in the overlap
psd=pwelch(x,hann(nfft),noverlap,nfft,fs); % PSD
df=1/(nfft*dt);                    % frequency resolution
f = 0:df:fs/2;                     % frequency vector

%% calculation of PSD
Na = 32;                            
nfft=round(N/Na);
noverlap=round(nfft/2);
psd2=pwelch(x,hann(nfft),noverlap,nfft,fs);
df2=1/(nfft*dt); 
f2 = 0:df2:fs/2; 

%% Parseval's theorem
MSV=sum((x).^2)/T*dt               % mean square value
STD=sqrt(MSV);                     % standard deviation
PSD=trapz(psd2)*df2                % area under the PSD plot   
amp=PSD/(fs/2)                     % average value of PSD 


%% Plot of the results
figure (1)                         % time domain plot
t1=0.001:0.1:T;
plot(t,x,'linewidth',2,'Color',[.6 .6 .6])
hold on
plot(t1,STD*t1./t1,'--k','linewidth',3)
set(gca,'fontsize',16)
axis([0,5,-4,4])
xlabel('time(s)');
ylabel('displacement (m)');
grid, axis square

figure (2)                              % frequency domain plot
fa=0:600;
plot(f,10*log10(psd),'linewidth',2,'Color',[.6 .6 .6])
hold on
plot(f,10*log10(psd),'o','linewidth',2,'Color',[.6 .6 .6])
hold on
plot(f2,10*log10(psd2),'k','linewidth',2)
hold on
plot(f2,10*log10(psd2),'sk','linewidth',2)
plot(fa,10*log10(fa./fa*amp),'--k','linewidth',3)
xticks([0 20 40 60 80 100])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('PSD (dB ref 1 m^2/Hz)');
axis([0,fs/2,-30,-15]),grid, axis square




