%% Exercise 2.5
clear all;

%% Parameters
m=1;                                % mass
k=1e4;                              % stiffness
z=0.01;                             % damping ratio
c=2*z*sqrt(k*m);                    % damping coefficient

%% Frequency vector
f=0:1:20;                           % frequency vector
w=2*pi*f;

%% Dynamic Stiffness
K=k-w.^2+j*w*c;                     % dynamic stiffness

%% Calculations
ff=f.^2;
p = polyfit(ff,real(K),1);
stiffness=p(2)
mass=-p(1)/(2*pi)^2
q=polyfit(f,imag(K),1);
damping=q(1)/(2*pi)/(2*sqrt(mass*stiffness))

%% Plot results
figure (1)
plot(ff,real(K),'ok','markersize',8,'linewidth',2)
set(gca,'fontsize',16)
xlabel('frequency^2 (Hz^2)');
ylabel('real(dynamic stiffness) (N/m)');
axis([0,400,-1e4,1e4])
grid;axis square

figure (2)
plot(f,imag(K),'ok','markersize',8,'linewidth',2)
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('imag(dynamic stiffness) (N/m)');
axis([0,20,0,300])
grid;axis square

