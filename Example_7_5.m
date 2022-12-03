%% virtual experiment on a vibration isolator II
clear all 
%% Parameters
% suspended mass
m = 1;                              % mass           
% isolator
k=10000;                            % stifness  
z = 0.1; c = 2*z*sqrt(m*k);         % damping 

%% calculation of the transmitted forces
%slow rate of change of freq.
n=1000;                             % number of cycles
[fmax,fe,ft,t]=calc(m,k,c,n); 
fe1=fe;ft1=ft;t1=t;

%medium rate of change of freq.     
n=200;                              % number of cycles
[fmax,fe,ft,t]=calc(m,k,c,n); 
fe2=fe;ft2=ft;t2=t;

%fast rate of change of freq.
n=20;                               % number of cycles
[fmax,fe,ft,t]=calc(m,k,c,n); 
fe3=fe;ft3=ft;t3=t;

%% plot the results
plot(t1,ft1,'linewidth',3,'Color',[.6 .6 .6])
hold on
plot(t2,ft2,'linewidth',3,'Color',[.4 .4 .4])
hold on
plot(t3,ft3,'linewidth',3,'Color',[.2 .2 .2]); grid
set(gca,'fontsize',24)
axis([0,12,1.1*min(ft1),1.1*max(ft1)])
axis([0,12,-15,15])
xlabel('time (s)');
ylabel('transmitted force (N)');

fn=1/(2*pi)*sqrt(k/m)
fmax

function [fmax,fe,ft,t]=calc(m,k,c,n)  % function to calculate the transmitted force
fs=2000;dt=1/fs;                      % sampling frequency/. time resolution
fmax=100;                             % operational frequency
T=n/fmax;t1=0:dt:T;                   % time vector  
a=2*pi*fmax/(2*T);                    % coefficient
ff1=sin(a*t1.^2);                     % chirp signal
t2=0:dt:15-T;                         % time vector
ff2=sin(n*pi+2*pi*fmax*t2);           % steady-state force
fe=[ff1, ff2];                        % total force
N=length(fe);
t=0:dt:(N-1)*dt; Tm=max(t);           % time vector

df=1/(N*dt);                          % frequency resolution
ff = 0:df:fs/2;                       % frequency vector
w=2*pi*ff;ft=ff;
%% calculation of the transmitted force
Tf=(k+j*w*c)./(k-w.^2*m+j*w*c);           % force transmissibility
TTf=[Tf,fliplr(conj(Tf(1:length(Tf)-1)))];
ttf=fs*ifft(TTf);                         % IDFT of the transmissibility
ft=conv(real(ttf),fe)/fs;                 % transmitted force   
ft=ft(1:length(fe));                      % transmitted force 
end



