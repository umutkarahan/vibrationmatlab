clear all

%% Time vectors
dt=0.05; dt1=0.001;              % time resolution in seconds
T=1;                             % duration of time signal
t=0:dt:T; t1=0:dt1:T;            % time vectors
f=1;                             % frequency in Hz
w=2*pi*f;                        % frequency in rad/s

%% velocity
v=sin(w*t1);                     % velocity – high time resolution 
vm=sin(w*t);                     % velocity – low time resolution

%% differentiation
a=w.*cos(w*t1);                  % acceleration - high time resolution
am=diff(vm)/dt;                  % numerical differentiation
tt=dt/2:dt:T-dt/2;               % define time vector (one point less)

%% integration
x=-1/w*cos(w*t1);                % displacement - high time resolution
xt=cumtrapz(vm)*dt + x(1);       % numerical integration

%% plot the results
figure                           % displacement as a function of time
plot(t,xt,'ok',t1,x,'k','linewidth',2, 'MarkerSize',10)
set(gca,'fontsize',16)
xticks([0 0.2 0.4 0.6 0.8 1]) 
axis([0 1 -0.2 0.2])
grid;axis square
xlabel('time (s)');
ylabel('displacement (m)');

figure                            % velocity as a function of time
plot(t1,v,'k',t,vm,'ok','linewidth',2, 'MarkerSize',10)
set(gca,'fontsize',16)
xticks([0 0.2 0.4 0.6 0.8 1]) 
axis([0 1 -1 1])
grid;axis square
xlabel('time (s)');
ylabel('velocity (m/s)');

figure                            % acceleration as a function of time
plot(t1,a,'k',tt,am,'ok','linewidth',2, 'MarkerSize',10)
set(gca,'fontsize',16)
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([-8 -6 -4 -2 0 2 4 6 8])
axis([0 1 -8 8])
grid;axis square
xlabel('time (s)');
ylabel('acceleration (m/s^2)');