%% program to calculate IRFs for MDOF systems
clear all

%% parameters
m1=1;m2=1;m3=1;                               % masses
k1=1e4;k2=1e4;k3=1e4;k4=0*1e4;                % stiffness
M=[m1 0 0; 0 m2 0; 0 0 m3];                   % mass matrix
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];  % stiffness matrix
C=5e-4*K;                                     % damping matrix

%% Modal matrices
[V W]= eig (inv(M)*K);                         % calculation of eigenvalues and eigenvectors
R=sqrt(W);                                     % Calculation of natural frequencies
Mm=V'*M*V;                                     % modal mass matrix
Km=V'*K*V;                                     % modal stiffness matrix
Cm=V'*C*V;                                     % modal damping matrix
zeta=Cm/(2*sqrt(Km*Mm));                       % modal damping ratios
wd=R*sqrt(1-zeta(1,1)^2);                      % damped natural frequencies

w1=R(1,1);w2=R(2,2);w3=R(3,3);                 % undamped natural frequencies
wd1=wd(1,1);wd2=wd(2,2);wd3=wd(3,3);           % damped natural frequencies
z1=zeta(1,1);z2=zeta(2,2);z3=zeta(3,3);        % modal damping ratios

%% forced vibration
fs=1000;df=0.001;dt=1/fs;                      % time and frequency parameters                                 
n=0;
for f=0:df:fs/2                             % frequency vector
    w=2*pi*f;
    n=n+1;
    An=inv(K-w.^2*M+j*w*C);                 % Matrix of FRFs
    An11(n)=An(1,1);                        % X1/F1
    An21(n)=An(2,1);                        % X2/F1
    An31(n)=An(3,1);                        % X3/F1

    % Modal responses
    Hm=inv(Km-w^2*Mm+j*w*Cm);
    % mode 1
    F=[1 0 0]';                             % force vector
    Q=[1 0 0; 0 0 0;0 0 0];                 % picks out first mode
    Xa=V*Q*Hm*V'*F;                         % calculates first mode in FRFs
    % mode 2
    Q=[0 0 0; 0 1 0;0 0 0];                 % picks out second mode
    Xb=V*Q*Hm*V'*F;                         % calculates second mode in FRFs
    % mode 3
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

    % X3/F1                                  % individual modal responses for X3/F1
    X31a(n)=Xa(3);
    X31b(n)=Xb(3);
    X31c(n)=Xc(3);
end

%% IRFs
HH=An11;HH1=X11a;HH2=X11b;HH3=X11c;              % FRFs for X1/F1
[h h1 h2 h3]=IRF(fs,HH,HH1,HH2,HH3);             % function to calculate IRFs
t=0:dt:(length(h)-1)*dt;                         % time vector
h11=h;h11a=h1;h11b=h2;h11c=h3;                   % overall and modal IRFs 

A1=1/wd1*exp(-z1*w1*t).*sin(wd1*t);              % part of mode 1 IRF
A2=1/wd2*exp(-z2*w2*t).*sin(wd2*t);              % part of mode 2 IRF
A3=1/wd3*exp(-z3*w3*t).*sin(wd3*t);              % part of mode 3 IRF

h11at=V(1,1)^2*A1;                               % mode 1 IRF for X1/F1
h11bt=V(1,2)^2*A2;                               % mode 2 IRF for X1/F1
h11ct=V(1,3)^2*A3;                               % mode 3 IRF for X1/F1
h11t=h11at+h11bt+h11ct;                          % overall IRF for X1/F1 

HH=An21;HH1=X21a;HH2=X21b;HH3=X21c;              % FRFs for X2/F1
[h h1 h2 h3]=IRF(fs,HH,HH1,HH2,HH3);             % function to calculate IRFs
h21=h;h21a=h1;h21b=h2;h21c=h3;                   % overall and modal IRFs 

h21at=V(1,1)*V(2,1)*A1;                          % mode 1 IRF for X2/F1
h21bt=V(1,2)*V(2,2)*A2;                          % mode 2 IRF for X2/F1
h21ct=V(1,3)*V(2,3)*A3;                          % mode 2 IRF for X2/F1
h21t=h21at+h21bt+h21ct;                          % overall IRF for X2/F1 


HH=An31;HH1=X31a;HH2=X31b;HH3=X31c;              % FRFs for X3/F1
[h h1 h2 h3]=IRF(fs,HH,HH1,HH2,HH3);             % function to calculate IRFs
h31=h;h31a=h1;h31b=h2;h31c=h3;                   % overall and modal IRFs 

h31at=V(1,1)*V(3,1)*A1;                          % mode 1 IRF for X3/F1
h31bt=V(1,2)*V(3,2)*A2;                          % mode 2 IRF for X3/F1
h31ct=V(1,3)*V(3,3)*A3;                          % mode 3 IRF for X3/F1
h31t=h31at+h31bt+h31ct;                          % overall IRF for X3/F1 

%% figures
figure
plot(t,h11,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h11t,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
ylabel('IRF (m/Ns)')

figure
plot(t,h11a,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h11at,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
ylabel('Mode 1 IRF (m/Ns)')

figure
plot(t,h11b,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h11bt,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
ylabel('Mode 2 IRF (m/Ns)')
%%
figure
plot(t,h11c,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h11ct,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
xlabel('time (s)'),ylabel('Mode 3 IRF (m/Ns)'),
%%
figure
plot(t,h21,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h21t,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)

figure
plot(t,h21a,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h21at,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
figure
plot(t,h21b,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h21bt,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
figure
plot(t,h21c,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h21ct,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
xlabel('time (s)')

figure
plot(t,h31,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h31t,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
figure
plot(t,h31a,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h31at,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
figure
plot(t,h31b,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h31bt,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
figure
plot(t,h31c,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(t,h31ct,'--k','linewidth',3)
axis square; grid; axis([0,2,-7e-3,9e-3])
set(gca,'fontsize',16)
xlabel('time (s)')


function [h h1 h2 h3]=IRF(fs,HH,HH1,HH2,HH3)           % function to calculate IRFs
    Hd=[HH fliplr(conj(HH))];
    H=Hd(1:length(Hd)-1);
    h=fs*ifft(H);
    
    H1d=[HH1 fliplr(conj(HH1))];
    H1=H1d(1:length(H1d)-1);
    h1=fs*ifft(H1);
    
    H2d=[HH2 fliplr(conj(HH2))];
    H2=H2d(1:length(H2d)-1);
    h2=fs*ifft(H2);
    
    H3d=[HH3 fliplr(conj(HH3))];
    H3=H3d(1:length(H3d)-1);
    h3=fs*ifft(H3);
end