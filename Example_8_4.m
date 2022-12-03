clear all

%% parameters
m1=1;m2=1;m3=1;                               % masses
k1=1e4;k2=1e4;k3=1e4;k4=0*5e3;                % stiffness
M=[m1 0 0; 0 m2 0; 0 0 m3];                   % mass matrix
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];  % stiffness matrix
C=1e-4*K;                                     % damping matrix

%% Modal matrices
[V W]= eig (inv(M)*K);                         % calculation of eigenvalues and eigenvectors
R=sqrt(W)/(2*pi);                              % Calculation of natural frequencies
Mm=V'*M*V;                                     % modal mass matrix
Km=V'*K*V;                                     % modal stiffness matrix
Cm=V'*C*V;                                     % modal damping matrix
zeta=Cm/(2*sqrt(Km*Mm))                        % modal damping ratios

%% forced vibration
n=0;
for f=0:0.001:50                               % frequency vector
    w=2*pi*f;
    n=n+1;
    An=inv(K-w.^2*M+j*w*C);
    An11(n)=An(1,1);
    An21(n)=An(2,1);
    An31(n)=An(3,1);

    % Modal responses
    Hm=inv(Km-w^2*Mm+j*w*Cm);
    % mode 1
    F=[1 0 0]';                             % force vector
    Q=[1 0 0; 0 0 0;0 0 0];                 % picks out first mode
    Xa=V*Q*Hm*V'*F;                         % calculates first mode in FRFs
    %mode 2
    Q=[0 0 0; 0 1 0;0 0 0];                 % picks out second mode
    Xb=V*Q*Hm*V'*F;                         % calculates second mode in FRFs
    %mode 3
    Q=[0 0 0; 0 0 0;0 0 1];                 % picks out third mode
    Xc=V*Q*Hm*V'*F;                         % calculates third mode in FRFs

    % X1/F1                                 % individual modal responses for X1/F1
    X11a(n)=Xa(1);
    X11b(n)=Xb(1);
    X11c(n)=Xc(1);

    % X2/F1                                 % individual modal responses for X2/F1
    X21a(n)=Xa(2);
    X21b(n)=Xb(2);
    X21c(n)=Xc(2);

    %X3/F1                                  % individual modal responses for X3/F1
    X31a(n)=Xa(3);
    X31b(n)=Xb(3);
    X31c(n)=Xc(3);

end

f=0:0.001:50;
figure
plot(f,20*log10(abs(An11)),'k','linewidth',4)
hold on
plot(f,20*log10(abs(X11a)),'--','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X11b)),':','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X11c)),'linewidth',3,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*(angle(X11a)),'--','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X11b)),':','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X11c)),'linewidth',4,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-200,200])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure
plot(f,20*log10(abs(An21)),'k','linewidth',4)
hold on
plot(f,20*log10(abs(X21a)),'--','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X21b)),':','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X21c)),'linewidth',3,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*(angle(X21a)),'--','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X21b)),':','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X21c)),'linewidth',4,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-200,200])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure
plot(f,20*log10(abs(An31)),'k','linewidth',4)
hold on
plot(f,20*log10(abs(X31a)),'--','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X31b)),':','linewidth',3,'color',[0.5 0.5 0.5])
hold on
plot(f,20*log10(abs(X31c)),'linewidth',3,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-130,-30])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('amplitude (dB ref 1 m/N)');

figure
plot(f,180/pi*(angle(X31a)),'--','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X31b)),':','linewidth',4,'color',[0.5 0.5 0.5])
hold on
plot(f,180/pi*(angle(X31c)),'linewidth',4,'color',[0.5 0.5 0.5])
axis square; grid; axis([0,40,-200,200])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');