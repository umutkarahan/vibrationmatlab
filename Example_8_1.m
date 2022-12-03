clear all

%% parameters
m1=1;m2=1;m3=1;                               % masses
k1=1e4;k2=1e4;k3=1e4;k4=0*5e3;                % stiffness
M=[m1 0 0; 0 m2 0; 0 0 m3];                   % mass matrix
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];  % stiffness matrix
C=1e-4*K;                                     % damping matrix

%% forced vibration
n=0;
for f=0:0.001:50                         % frequency vector
 w=2*pi*f;
 n=n+1;
 An=inv(K-w.^2*M+j*w*C);                 % calculate matrix of FRFs for each frequency
 An11(n)=An(1,1);                        % point receptance
 An21(n)=An(2,1);                        % transfer receptance
 An31(n)=An(3,1);                        % transfer receptance
end

f=0:0.001:50;
figure
plot(f,20*log10(abs(An11)),'k','linewidth',4)
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*unwrap(angle(An11)),'k','linewidth',4)
axis square; grid; axis([0,40,-600,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure
plot(f,20*log10(abs(An21)),'k','linewidth',4)
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*unwrap(angle(An21)),'k','linewidth',4)
axis square; grid; axis([0,40,-600,0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure
plot(f,20*log10(abs(An31)),'k','linewidth',4)
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*unwrap(angle(An31)),'k','linewidth',4)
axis square; grid; axis([0,40,-600,0])
ylabel('phase (degrees)');
set(gca,'fontsize',16)
xlabel('frequency (Hz)');