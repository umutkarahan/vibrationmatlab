clear all 

%% parameters
m = 1;                              % [kg]           see MATLAB Example 3.1
k=10000;                            % [N/m]  
z = 0.1; c = 2*z*sqrt(m*k);         % [Ns/m] 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % [rad/s]

%% time and frequency parameters
T=0.5;                              % [s]          T needs to be changed to 1s to remove wraparound 
fs=500; dt=1/fs; t=0:dt:T;             

%% IRF
ho=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);             % IRF%
%ho((length(t)-1)/2:length(t))=0;                % this needs to be uncommented when T=1

%% force
N=1;N=200;
N1=25; N2=length(t)-N1-N;
fo=[zeros(1,N1) ones(1,N) zeros(1,N2)];          % time hisory of the force

%% convolution
y = dt*conv(ho,fo);                              % response by convolution
tt=0:dt:2*T;                                     % time vector for the response

%% circular convolution
H=fft(ho)*dt;                                    % FRF
F=fft(fo)*dt;                                    % force spectrum
YC=H.*F;                                         % response spectrum
yc=ifft(YC)*fs;                                  % response by circular convolution

%% plot the results
figure
plot(t,fo,'k','linewidth',3)                     % force
ylim([0,2])
set(gca,'fontsize',16)
axis square;grid;
axis([0,T,0,2])
xticks([0 .2 .4 .6 .8 1])
xlabel('time (s)');
ylabel('force (N)');

figure
plot(t,ho,'k','linewidth',3)                      % displacement IRF
set(gca,'fontsize',16)
axis square;grid;
axis([0,T,-1e-2, 1e-2])
%xticks([0 .1 .2 .3 .4 .5])
xticks([0 .2 .4 .6 .8 1])
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');

figure
gr=[.6 .6 .6];                                     
plot(tt,y,'k','linewidth',3);grid on; hold on;     % response by convolution 
plot(t,yc,'--k','linewidth',3,'color',gr);         % response by circular convolution
set(gca,'fontsize',16)
axis square;
xlabel('time (s)');
ylabel('displacement (m)');
