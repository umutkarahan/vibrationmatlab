clear all 

%% Parameters
m = 1;                              % see Matlab example 3.1           
k=10000;                              
z = 0.1; c = 2*z*sqrt(m*k);          
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    

%% Time and frequency parameters
Th=1;                                % see Matlab example 3.1
fs=400;                              
dt=1/fs; t=0:dt:Th;                  
df=1/Th; f=0:df:fs;                  
Tw=.2; tw=0:dt:Tw;                   
dfw=1/Tw; fw=0:dfw:fs;               
%% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);           % IRF long window
ht=1/(m*wd)*exp(-z*wn*tw).*sin(wd*tw);        % IRF short window
w=[ones(1,length(tw)) zeros(1,(length(t)-(length(tw))))];  % long window
hw=h.*w;                                      % truncated FRF

%% plot IRFs
figure (1)
plot(t,h,'k','linewidth',4)                   % IRF long window
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');
grid;axis square
figure (2)                                    % IRF long window truncated IRF
plot(t,hw,'k','linewidth',4)
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');
grid;axis square
figure (3)
plot(t,h,'k','linewidth',4)                  % IRF truncated short window
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');
grid;axis square
axis([0 0.2 -inf inf])

%% Calculation of DFT
HW=dt*fft(hw);                               % DFT – no truncation
H=dt*fft(h);                                 % DFT – truncation, short window
HTT=dt*fft(ht);                              % DFT – truncation, long window

%% Theoretical FRF
dff=0.1; fr=0:dff:fs/2;             % [Hz]
w=2*pi*fr;                          % [rad/s]
HH=1./(k-w.^2*m+j*w*c);             % [m/N]

%% Plot results
figure (4)                                                  % modulus
plot(f,20*log10(abs(HW)),'linewidth',3,'Color',[.4 .4 .4])
hold on
plot(fr,20*log10(abs(HH)),'k','linewidth',4)
hold on
plot(f,20*log10(abs(H)),'--k','linewidth',4,'Color',[.6 .6 .6])
hold on
plot(fw,20*log10(abs(HTT)),'k','linewidth',2)
axis([0 fs/2 -130 -60])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
grid;axis square

figure (5)                                                  % phase
plot(f,180/pi*unwrap(angle(HW)),'linewidth',3,'Color',[.4 .4 .4])
hold on
plot(fr,180/pi*unwrap(angle(HH)),'k','linewidth',4)
hold on
plot(f,180/pi*unwrap(angle(H)),'--k','linewidth',4,'Color',[.6 .6 .6])
hold on
plot(fw,180/pi*unwrap(angle(HTT)),'k','linewidth',2)
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('phase (^o)');
grid;axis square
axis([0 fs/2 -250 0])
