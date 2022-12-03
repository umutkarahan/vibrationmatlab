
clear all
%% time and frequency data
dt=0.0001;                  % time resolution            
T = 3;                      % time duration
t = 0:dt:T;                 % time vector
f=10;                       % frequency
w=2*pi*f;

%% signals and the envelope
x=1*sin(w*t);               % signal x
y=0.9*sin(0.95*w*t);        % signal y
z=x+y;                      % sum of signals x and y
zh=hilbert(z);              % Hilbert transform of z

%% plot the results
plot(t,z,'k','linewidth',2)         % plot z and its envelope
hold on
plot(t,abs(zh),'--k',t,-abs(zh),'--k','linewidth',2) 
axis([0 3,-3,3])
axis square; grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement (m)');


