%% Exercise 2.3
clear all

%% Parameters
m=1;                                % mass
k=1e4;                              % stiffness
wn=sqrt(k/m);                       % natural frequency
z=0.01;                             % damping ratio
c=2*z*wn*m;                         % damping coefficient

%% Frequency vector
df=0.001;                           % frequency resolution
F=100;                              % maximum frequency
f=0:df:F;                           % frequency vector
w=2*pi*f;                           % circular frequency

%% Receptance FRF
H=1./(k-w.^2+j*w*c);                % receptance FRF

%% Plot results
figure (1)                                                    % modulus
semilogx(f,20*log10(abs(H)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,20*log10(1./w.^2),'--k','linewidth',2)
hold on
semilogx(f,20*log10(f./f*1/k),':k','linewidth',2)
hold on
semilogx(wn/(2*pi),20*log10(1/(wn*c)),'xk','linewidth',6)
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
grid;axis square
axis([1e0,100,-110,-40])

figure (2)                                                    % phase
semilogx(f,180/pi*angle(H),'k','linewidth',4)
grid;axis square
axis([1e0,100,-200,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase angle (degrees)');

figure (3)                                                    % real part
plot(f,real(H),'k','linewidth',4)
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('real(receptance) (m/N)');
grid;axis square
axis([10,20,-4e-3,4e-3])

figure (4)                                                    % imaginary
plot(f,imag(H),'k','linewidth',4)
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('imag(receptance) (m/N)');
grid;axis square
axis([10,20,-6e-3,0])

figure (5)                                                    % Nyquist
plot(real(H),imag(H),'k','linewidth',4)
set(gca,'fontsize',16)
xlabel('real (receptance) (m/N)');
ylabel('imag(receptance) (m/N)');
grid;axis square
axis([-3e-3,3e-3,-6e-3,0])
