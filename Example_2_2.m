%% Exercise 2.2
clear all 
 
%% Parameters
m=1;                                      % mass
k=10000;                                  % stiffness
wn=sqrt(k/m);                             % undamped natural frequency
z=0.01;                                   % damping ratio
c=2*z*wn*m;                               % damping coefficient
wd=sqrt(1-z^2)*wn;                        % damped natural frequency
Td=2*pi/wd;                               % damped natural period
%% Time Vector
dt=0.001;                                 % time resolution
T=5;                                      % duration of time signal
t=0:dt:T;                                 % time vector

%% Normalised displacement IRF
h=exp(-z*wn*t).*sin(wd*t);                % normalised IRF

%% Calculations
a=hilbert(h);                             % create analytic signal
env=log(abs(a));                          % calculate log of the envelope
t1=t(500:3500);env1=env(500:3500);        % set a specific time range
p=polyfit(t1,env1,1);                     % fit a straight line in the range t1
grad=-p(1);                               % calculate the gradient of the line
gamma=2*pi/(grad*Td);                     % calculate the constant
Est_z=1/(1+gamma)                         % damping ratio estimate 


%% Plot results
figure(1);
plot(t,h,'linewidth',2,'Color',[1 1 1]*0.6);
hold on;
plot(t,abs(a),'k',t,-abs(a),'k','linewidth',2)
grid;axis square
set(gca,'fontsize',16)
xlim([0,5])
xlabel('time (s)');
ylabel('normalised displacement');

figure(2);
plot(t,env,'k','linewidth',2);
hold on
plot(t1,env1,'k','linewidth',8,'Color',[.6 .6 .6])
grid;axis square
axis([0,5,-5,0])
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('ln(envelope)');


