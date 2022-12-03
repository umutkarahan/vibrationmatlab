clear all

%% Parameters
m = 1;                              % mass           
k = 10000;                          % stiffness  
z = 0.08;                           % damping ratio
c = 2*z*sqrt(m*k);                  % damping coefficient 
wn=sqrt(k/m);                       % undamped natural frequency
wd=sqrt(1-z^2)*wn;                  % damped natural frequency

%% Time and frequency parameters
T=2;                                % duration of time signal
fs=400;                             % sampling frequency
dt=1/fs;                            % time resolution
t=0:dt:T;                           % time vector
df=1/T;                             % frequency resolution
f=0:df:fs;                          % frequency vector
%% Theoretical IRF
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % IRF

%% Calculation of DFT
H=dt*fft(h);                        % calculation of receptance FRF

%% Theoretical FRF
dff=df;                             % frequency resolution
fr=0:dff:fs/2;                      % frequency vector
w=2*pi*fr;                          % frequency in rad/s
HH=1./(k-w.^2*m+j*w*c);             % theoretical receptance FRF

%% IFT
Hd=[HH fliplr(conj(HH))];
Hdd=Hd(1:length(Hd)-1);    
hd=fs*ifft(Hdd);

%% Plot results
figure (1)
plot(f,20*log10(abs(Hdd)),'linewidth',4,'Color',[0.6 0.6 0.6]) % modulus
hold on
plot(f,20*log10(abs(H)),'--k','linewidth',2)
grid;axis square
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
%%
figure (2)
plot(f,180/pi*angle(Hdd),'linewidth',4,'Color',[0.6 0.6 0.6]) % phase
hold on
plot(f,180/pi*angle(H),'--k','linewidth',2)
grid;axis square
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
%%
figure (3)
plot(t,hd,'linewidth',4,'Color',[0.6 0.6 0.6])                % IRF
hold on
plot(t,h,'--k','linewidth',2)
axis([0,1,-0.01,0.01])
grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('IRF (m/Ns)');




 