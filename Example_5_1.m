clear all 

%% time and frequency parameters
fn=2; fs=40;                           % excitation and sampling frequencies
Nc=2.4; A=1;                           % number of cycles; amplitude of excitation 
T=Nc/fn; dt=1/fs; t=0:dt:T-1/fs;       % time vector
N=length(t);                           % number of points 
TT=N*dt;

%% sine wave
x = A*sin(2*pi*fn*t);                  % sine wave

%% calculation of PSD
xdft = fft(x)*dt;                      % DFT of sine wave 
psd = abs(xdft).^2/TT;                 % PSD of sine wave
df=1/TT;                               % frequency resolution
f=0:df:(N-1)*df;                       % frequency vector
 
%% Parseval's theorem
MSV=sum((x).^2)/(TT)*dt                      % mean square value
PSD=2*trapz(psd(1:round(N/2+1)))*df          % area under the PSD plot

%% Plot of the results
figure (1)
plot(t,x,'k','linewidth',2)
hold on
plot(t,x,'ok','linewidth',2,'MarkerSize',10)
axis([0,inf,-1,1])
xticks([0 .2 .4 .6 .8 1 1.2])
set(gca,'fontsize',16)
xlabel('time(s)');
ylabel('displacement (m)');
grid;axis square

figure (2)
plot(f,psd,'k','linewidth',2)
hold on
plot(f,psd,'ok','linewidth',2,'MarkerSize',10)
axis([0 fs, 0, 0.25])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('PSD (m^2/Hz)');
grid;axis square



