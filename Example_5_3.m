clear all

%% time and frequency parameters
fs=40;                             % sampling frequency 
T = 8;dt=1/fs;t=0:dt:T;            % chirp duration; time resolution; time vector                  
N=length(t);A=2;                   % number of points; amplitude of chirp
f1=1; f2=5;                        % upper and lower frequencies
TT=N*dt;                            

%% chirp
a=2*pi*(f2-f1)/(2*T); b=2*pi*f1;   % coefficients
x=A*sin(a*t.^2+b*t);               % chirp signal

%% calculation of ESD
X = fft(x)/fs;                        % DFT of chirp
esd = abs(X).^2;                      % ESD of chirp
df=1/(N*dt);                          % frequency resolution
f = 0:df:fs/2 ;                       % frequency vector
esds=2*abs(X(1:round(N/2))).^2;       % ESD of chirp (single-sided)
esds(1)=abs(X(1)).^2;                 

%% Parsevals theorem
e=trapz((x).^2)*dt                    % energy in the time domain
E=trapz(esds)*df                      % energy in the frequency domain
amp=A^2/2*T/(f2-f1)                   % average amplitude 

%% Plot the results
figure(1)                             % time domain
plot(t,x,'k','linewidth',2)
hold on
plot(t,x,'ok','linewidth',2,'MarkerSize',10)
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time(s)');
ylabel('displacement (m)');
grid

figure(2)                             % frequency domain
plot(f,esds,'k','linewidth',2)
hold on
plot(f,esds,'ok','linewidth',2,'MarkerSize',10)
hold on
plot(f,amp*f./f,'--k','linewidth',2)
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('ESD (m/Hz)^2');
axis([0,6,0,inf]),grid






