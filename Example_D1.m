clear all 

%% Variables

m = 1;                  % mass
k = 1000;               % stiffness
z = 0.05;               % damping ratio
c = 2*z*sqrt(m*k);      % damping factor
A = [0 1; -k/m -c/m];   % system matrices
B = [0; 1/m];           %                       

dt = 0.005;             % time resolution in seconds
fs = 1/dt;              % sampling frequency
T = 10;                 % duration of time signal
df=1/T;                 % frequency resolution

t = 0:dt:T;             % time vector
f = 0:df:fs;            % frequency vector
L = length(t);          % number of samples

%% definition of an input force time series
f1 = zeros(size(t));
Tn = T*5/12;                                                         % time duration of the force pulse 
Ln = round(Tn/dt);                                                   % number of samples in the force pulse
f = [zeros(1,100), ones(1,Ln), zeros(1,(length(t)-(Ln+100)))];       % force time history

%% Runge-Kutta Method
n=t;                                                                 %dummy variable for the look up table used in the function interp1
[t,x] = ode45(@(t,x) pulse(t,x,A,B,f,n),t,[0 0]);                    % response due to the force pulse 

%% plot the results
figure
plot(t,f,'k','linewidth',2);
axis([0,10,0,1.2])
xticks([0 2 4 6 8 10])
axis square; grid
set(gca,'fontsize',16)
xlabel('time (s)')
ylabel('excitation force (N)')
figure
plot(t,x(:,1),'k','linewidth',2)
xticks([0 2 4 6 8 10])
axis square; grid
set(gca,'fontsize',16)
xlabel('time (s)')
ylabel('displacement (m)')
figure
plot(t,x(:,2),'k','linewidth',2)
xticks([0 2 4 6 8 10])
axis square; grid
set(gca,'fontsize',16)
xlabel('time (s)')
ylabel('velocity (m/s)')

%% function
function dx = pulse(t,x,A,B,f,n)
    force = interp1(n,f,t);
    dx=A*x+B*force;
end