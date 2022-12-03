%% calculation of the output from a force impulse
clear all 
%% Parameters
m = 1;                              % mass           
k=1000;                             % stifness  
z = 0.1; c = 2*z*sqrt(m*k);         % damping 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % natural frequency
fn=wn/(2*pi);                       % natural frequency
Tn=1/fn;                            % natural period 

%% Example 6.1a Impulsive force
TT=2;                               % time duratiom
fs=100;                             % samping frequency
dt=1/fs; t=0:dt:TT-dt;              % time vector
f = zeros(1,length(t));f(12) = 1;   % force vector 
[xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs); 
plots(t,f,TT,xc,xf,xn)
%% Example 6.1b. Half-sine pulse
Tc=Tn*0.5;                          % [s]
fs=100;                             % [Hz]
dt=1/fs; t=0:dt:Tc;                 % [s]
f = sin(pi*t/Tc);                   % half sine pulse
Nz=20*length(t);                    % Number of zeros to add
f = [f zeros(Nz,1)'];               % zero-padded force signal
t = 0:dt:(length(f)-1)*dt;          % time vector for extended signal
N=length(t);TT=N*dt;
[xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs); 
plots(t,f,TT,xc,xf,xn)
%% Example 6.1c. Chirp
% slow chirp
fs=1000;                           % sampling frequency 
T=60;dt=1/fs;t=0:dt:T;             % chirp duration; time resolution; time vector                  
f1=1;f2=10;                        % upper and lower frequencies

[f,t,TT]=chrp(f1,f2,t,T,dt);
[xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs); 
plots(t,f,TT,xc,xf,xn)

% fast chirp
T=1.25;t=0:dt:T;                   % chirp duration; time resolution; time vector                  
f1=1;f2=100;                       % upper and lower frequencies
[f,t,TT]=chrp(f1,f2,t,T,dt); 
[xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs); 
plots(t,f,TT,xc,xf,xn)

%% Example 6.1d. Random excitation
fs=100;                            % sampling frequency 
T=10;dt=1/fs;t=0:dt:T;             % signal duration; time resolution; time vector                  
fc = randn(1,length(t));           % random signal 
fc=fc-mean(fc);                    % set the mean to zero

f=[fc zeros(1,length(fc))];        % zero padded force signal
t=0:dt:2*T+dt;                     % time vector for extended signal
N=length(t);TT=N*dt;  
[xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs); 
plots(t,f,TT,xc,xf,xn)

%% Functions

function [xc,xf,xn]=calculate(f,m,k,c,wd,wn,z,t,fs)   % function to calculate the response
%% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);                   % IRF 

%% Convolution
xc = conv(h,f)/fs;                                    % displacement response   
xc = xc(1:length(f));                                 % displacement response   

%% Frequency domain method
H = fft(h)/fs;                                        % FRF
F = fft(f)/fs;                                        % DFT of force
X = H.*F;                                             % DFT of displacement response
xf = ifft(X)*fs;                                      % displacement response

%% Numerical integration
A=[0 1;-k/m -c/m];                                    % system matrix
B=[0; 1/m];                                           % system matrix
n=t;                                                  % dummy variable for the look up table used in the function interp1
[t,xn]=ode45(@(t,x) imp(A,B,x,f,t,n),t,[0 0]);        % response due to the impulse 
xn=xn(:,1);
end

function [f,t,TT] = chrp(f1,f2,t,T,dt)                % function for chirp          
a=2*pi*(f2-f1)/(2*T); b=2*pi*f1;                      % coefficients
fc=sin(a*t.^2+b*t);                                   % chirp signal
f=[fc zeros(1,length(fc))];                           % zero padded force signal
t=0:dt:2*T+dt;                                        % time vector for extended signal
N=length(t);TT=N*dt; 
end
    
function dxdt=imp(A,B,x,f,t,n)                        % used in numerical solution
f = interp1(n,f,t);
dxdt=A*x+B*f;
end

function plots(t,f,TT,xc,xf,xn)                       % function for plots
figure
subplot(2,1,1)                                        % force input
plot(t,f,'-k','linewidth',4),grid
set(gca,'fontsize',24)
axis([0,TT,1.1*min(f),1.1*max(f)])
xlabel('time (s)');
ylabel('force (N)');
subplot(2,1,2)                                        % displacement response
plot(t,xc,'linewidth',4,'Color',[.6 .6 .6]),grid
hold on
plot(t,xf,'--k','linewidth',4)
hold on
plot(t,xn,'k','linewidth',1)
set(gca,'fontsize',24)
axis([0,TT,1.1*min(xc),1.1*max(xc)])
xlabel('time (s)');
ylabel('displacement (m)');
 end


