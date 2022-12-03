%% Exercise 2.1
clear all 

%% Parameters
m=1;                                % mass
k=10000;                            % stiffness
wn=sqrt(k/m);                       % undamped natural frequency
z=0.1;                              % damping ratio (set to 0.001, 0.01, 0.1, 0.999)
c=2*z*wn*m;                         % damping coefficient
wd=sqrt(1-z^2)*wn;                  % damped natural frequency

%% Time Vector
dt=0.001;                           % time resolution
T=100;                              % duration of time signal
t=0:dt:T;                           % time vector

%% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % displacement IRF
 
%% Plot results
plot(t,h,'linewidth',2,'Color',[1 1 1]*0.2);    
axis([0,0.8,-0.01,0.01]);
grid;axis square
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');

