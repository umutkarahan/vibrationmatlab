clear all;

%% time and frequency parameters
fs=2000;                       % sampling frequency; excitation frequency 
T=0.01;                        % contact time
dt=1/fs;t=0:dt:T-dt;           % Time duration;time resolution; time vector
N=length(t);                   % Number of points
TT=N*dt;
 
%% half-sine pulse
x = sin(pi*t/T);               % half sine pulse

%% extended half-sine pulse signal
Nz=5*N;                            % Number of zeros to add
xe = [x zeros(Nz,1)'];             % extended signal
te = 0:dt:(length(xe)-1)*dt;       % time vector for extended signal
Ne=length(te);                     % Number of points in extended signal
TTe=Ne*dt;

%% Calculation of ESDs
X=fft(x)*dt;                       % DFT of half sine pulse
esd=conj(X).*X;                    % ESD of half sine pulse 
df=1/TT;fx=0:df:fs/2;              % frequency resolution; frequency vector
esds=2*abs(X(1:round(N/2)+1)).^2;  % ESD of half-sine pulse (single-sided)
esds(1)=abs(X(1)).^2;     

Xe=fft(xe)*dt;                          % DFT of extended signal
esde=conj(Xe).*Xe;                      % ESD of extended signal    
dfe=1/TTe;fxe=0:dfe:fs/2  ;             % frequency resolution; frequency vector
esdes=2*abs(Xe(1:round(Ne/2)+1)).^2;    % ESD of half-sine pulse (single-sided)
esdes(1)=abs(Xe(1)).^2;     

%% Theory
ft=0:1:250;w=2*pi*ft                     % frequency vector
Xt=T*2*pi*cos(w*T/2)./(pi^2-(w*T).^2);   % spectrum
esdts=2*abs(Xt).^2;                      % ESD of half-sine pulse (single-sided)
esdts(1)=abs(Xt(1)).^2;  

%% Parseval's theorem
e=trapz(x.^2)*dt                         % energy in the time domain
E1=2*trapz(esd(1:round(N/2+1)))*df       % energy in the frequency domain
E2=2*trapz(esde(1:round(Ne/2+1)))*dfe
%% plot the results
figure (1)
plot(te,xe,'k','linewidth',2)
hold on
plot(te,xe,'ok','linewidth',2,'Markersize',10)
hold on
plot(t,x,'ok','linewidth',2,'Markersize',10,'MarkerFaceColor',[.1 .1 .1])
grid;xlim([0 max(te)+dt])
xlabel('time(s)');
ylabel('displacement (m)');
set(gca,'fontname', 'arial','fontsize',24)

figure (2)
plot(ft,10*log10(abs(esdts)),'k','linewidth',2)
hold on
plot(fx,10*log10(abs(esds)),'ok','linewidth',2,'Markersize',10,'MarkerFaceColor',[.1 .1 .1])
hold on
plot(fxe,10*log10(abs(esdes)),'ok','linewidth',2,'Markersize',10)
grid,axis([0 250 -90, -40])
xlabel('frequency (Hz)');
ylabel('ESD (dB ref 1 m^2s/Hz)');
set(gca,'fontname', 'arial','fontsize',24)

