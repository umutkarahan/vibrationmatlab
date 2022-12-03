clear all

%% Parameters
m = 1;                              % mass           
k = 10000;                          % stiffness  
z = 0.001;                          % damping ratio
c = 2*z*sqrt(m*k);                  % damping coefficient 
wn=sqrt(k/m);                       % undamped natural frequency
wd=sqrt(1-z^2)*wn;                  % damped natural frequency

%% Time and frequency parameters
T=100;                              % duration of time signal
fs=1000    ;                        % sampling frequency
dt=1/fs;                            % time resolution
t=0:dt:T;                           % time vector
df=1/T;                             % frequency resolution
f=0:df:fs;                          % frequency vector
N=fs*T;                             % number of points - 1 
%% Theoretical IRF
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % IRF

%% Calculation of DFT
H=dt*fft(h);                        % calculation of receptance FRF

%% Theoretical FRF
dff=0.001;                          % frequency resolution
fr=0:dff:fs/2;                      % frequency vector
w=2*pi*fr;                          % frequency in rad/s
HH=1./(k-w.^2*m+j*w*c);             % theoretical receptance FRF

%% Plot results
figure (1)
plot(t,h,'linewidth',2,'Color',[0.6 0.6 0.6])                     % IRF
xticks([0 20 40 60 80 100])
grid;axis square
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');
%%
figure (2)
plot(fr,20*log10(abs(HH)),'linewidth',4,'Color',[0.6 0.6 0.6])    % modulus
hold on
plot(f,20*log10(abs(H)),'--k','linewidth',2)
xticks([0 200 400 600 800 1000])
grid;axis square
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
%%
figure (3)
plot(fr,180/pi*angle(HH),'linewidth',4,'Color',[0.6 0.6 0.6])     % phase
hold on
plot(f,180/pi*angle(H),'--k','linewidth',2)
xticks([0 200 400 600 800 1000])
grid;axis square
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
%%
figure (4)
plot(real(HH),imag(HH),'linewidth',4,'Color',[0.6 0.6 0.6])      % Nyquist
hold on
plot(real(H(1:N/2)),imag(H(1:N/2)),'--k','linewidth',2)
xticks([-0.03 -0.02 -0.01 0 0.01 0.02 0.03])
grid;axis square
axis([-0.03,0.03,-0.06,0])
set(gca,'fontsize',16)
xlabel('real\{receptance\} (m/N) ');
ylabel('imag\{receptance\} (m/N)');
%%
figure (5)
semilogx(fr,20*log10(abs(HH)),'linewidth',4,'Color',[0.6 0.6 0.6])  % modulus
hold on
semilogx(f(1:N/2),20*log10(abs(H(1:N/2))),'--k','linewidth',2)
xticks([1e0 1e1 1e2])
grid;axis square
axis([1,500,-140,-20])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
%%
figure (6)
semilogx(fr,180/pi*angle(HH),'linewidth',4,'Color',[0.6 0.6 0.6])   % phase
hold on
semilogx(f(1:N/2),180/pi*angle(H(1:N/2)),'--k','linewidth',2)
xticks([1e0 1e1 1e2])
grid;axis square
axis([1,500,-200,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');




 