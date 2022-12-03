%% Exercise 2.4
clear all;

%% Parameters
m=1;                                % mass
k=1e4;                              % stiffness
wn=sqrt(k/m);                       % natural frequency
z=0.01;                             % damping ratio
c=2*z*wn*m;                         % damping coefficient

%% Frequency vector
f=15:0.05:17;                       % frequency vector
w=2*pi*f;

%% Receptance FRF
H=1./(k-w.^2+j*w*c);                % receptance FRF

%% Plot results
figure (1)
plot(f,20*log10(abs(H)),'ok','markersize',8,'linewidth',2)   % modulus
axis([15,17,-60,-45])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
grid;axis square

figure (2)
plot(f,180/pi*angle(H),'ok','markersize',8,'linewidth',2)    % phase
grid;axis square
axis([-inf,inf,-180,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase angle (degrees)');

figure (3)
plot(real(H),imag(H),'ok','markersize',8,'linewidth',2)      % Nyquist
grid;axis square
axis([-3e-3,3e-3,-6e-3,0])
set(gca,'fontsize',16)
xlabel('real (receptance) (m/N)');
ylabel('imag(receptance) (m/N)');

