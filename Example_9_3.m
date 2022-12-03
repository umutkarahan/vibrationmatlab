%% FRF of a beam including shaker effects
clear all

%% parameters
% cantilver beam
E=69e9;                                       % Youngs modulus of aluminium (N/m^2)
rho=2700;                                     % density (kg/m^3);
l=0.75;b=0.02;d=0.01;S=b*d;I=b*d^3/12;        % geometrical parameters 
z=0.01;n=2*z;                                 % damping ratio and loss factor
Ed=E*(1+j*n);                                 % complex Young's modulus
m=rho*S*l;                                    % mass of the rod% Youngs modulus of aluminium (N/m^2)

%% frequency parameters
fs=3000;df=0.01;dt=1/fs;                      % frequency parameters
f=0.001:df:fs/2;                              % frequency vector
w=2*pi*f;

% beam FRFs
n=0;
for x=0.1:0.05:l;
n=n+1;
[ht,Htt] = calcFRF(E,I,rho,S,z,m,l,x,w,fs);   % calculate FRFs and IRFs
H(n,:)=Htt;
h(n,:)=ht;
end

%% shaker
ms = 0.1;                                     % mass  
ws=2*pi*10;                                   % natural frequency (rad/s)
ks=ws^2*ms;                                   % stiffness
zs=0.1; cs=2*zs*sqrt(ms*ks);                  % damping
Ks=ks-w.^2*ms+j*w*cs;                         % dynamic stiffness of shaker
%% input and output
T=120;t=0:dt:T;                               % signal duration; time vector                  
fe = randn(1,length(t));                      % random signal 
fe=fe-mean(fe);                               % set the mean to zero
K=1./H(1,:);                                  % dynamic stiffness of structure
Fe=K./(K+Ks);                                 % amount of force applied to structure
G=[Fe,fliplr(conj(Fe(1:length(Fe)-1)))];      % double sided spectrum
g=fs*ifft(G);                                 % IRF to determine force applied to structure
fes = conv(real(g),fe)/fs;                    % force applied to structure   
fes = fes(1:length(fe));                      % force applied to structure 
N=length(fes);
[dis] = calcResp(fes,fe,fs,h);                % calculation of the displacement responses

%% frequency domain calculations
[Sffe,Sff,Sww,He,CoHe,ff] = DisFRF(fes,fe,dis,N,fs,dt);      % calculation of freq. domain quantitities

%% plot the results
figure                                                      % force time history
plot(t,fes,'linewidth',2,'color',[0.6 0.6 0.6]) 
hold on
plot(t,fe,'linewidth',2,'color',[0.4 0.4 0.4]) 
axis square,axis([0,120,-8,8]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('force (N)');

p=14;                                                      % measurement position
figure
plot(t,dis(p,:),'linewidth',2,'color',[0.6 0.6 0.6])       % disp. time history  
axis square,axis([0,120,-30e-5,30e-5]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement (m)');

figure                                                      % force PSD
semilogx(ff,10*log10(abs(Sff)),'linewidth',4,'color',[0.4 0.4 0.4])
hold on
semilogx(ff,10*log10(abs(Sffe)),'linewidth',4,'color',[0.6 0.6 0.6])
axis square,axis([1,1000,-80,0]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('force PSD (dB ref 1 N^2/Hz)');
                                                    
figure                                                      % disp. time history
semilogx(ff,10*log10(abs(Sww(p,:))),'linewidth',4,'color',[0.5 0.5 0.5]) 
axis square,axis([1,1000,-200,-60]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('displ. PSD (dB ref 1 m^2/Hz)');

figure                                                     % FRF modulus
semilogx(ff,20*log10(abs(He(p,:))),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,20*log10(abs(H(p,:))),'--k','linewidth',4)
hold on 
axis square,axis([1,1000,-160,-40]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');

figure                                                     % FRF phase
semilogx(ff,180/pi*unwrap(angle(He(p,:))),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(f,180/pi*unwrap(angle(H(p,:))),'--k','linewidth',4)
axis square,axis([1,1000,-1000,0]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                                     % coherence
semilogx(ff,CoHe(p,:),'linewidth',4,'color',[0.5 0.5 0.5])
axis square,axis([1,1000,0,1]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('coherence');
%%
function [ht,Htt] = calcFRF(E,I,rho,S,z,m,l,x,w,fs);    % function to calculate FRF and IRF 
nmax=10;                                                % number of modes
kl(1)=1.87510;kl(2)=4.69409;kl(3)=7.85476;              % kl values 1-3
kl(4)=10.9956;kl(5)=14.1372;                            % kl values 4,5
n=6:nmax;                                      
kl(n)=(2*n-1)*pi/2;                                     % kl values > 5
for n=1:nmax;
A=(sinh(kl(n))-sin(kl(n)))./(cosh(kl(n))+cos(kl(n)));
xf=0.1;                                                 % force position
phi1=cosh(kl(n)*xf/l)-cos(kl(n)*xf/l)-...
A.*(sinh(kl(n)*xf/l)-sin(kl(n)*xf/l));

phi2=cosh(kl(n)*x/l)-cos(kl(n)*x/l)-...                 % response position
A.*(sinh(kl(n)*x/l)-sin(kl(n)*x/l));
wn=sqrt((E*I)./(rho*S))*(kl(n)).^2;                     % natural frequency 
 
Ht(n,:)=phi1*phi2./(m*(wn^2-w.^2+j*2*w*wn*z));          % FRF of each mode
end
Htt=sum(Ht);                                            % overall receptance FRF

%IRF
Hd=[Htt fliplr(conj(Htt))];                             % form the double-sided spectrum  
Hdt=Hd(1:length(Hd)-1);                                 % set the length of the FRF
ht=fs*ifft(Hdt);                                        % calculation of the IRF
end

function [dis] = calcResp(fes,fe,fs,h);                 % function to calculate displacement responses
for n=1:14;
dd = conv(real(h(n,:)),fes)/fs;                         % displacement response   
dis(n,:) = dd(1:length(fe)); 
end
end

function [Sffe,Sff,Sww,He,CoHe,ff] = DisFRF(fes,fe,dis,N,fs,dt)     % function to calculate frequency domain quantities
Na = 8;                                                             % number of averages
nfft=round(N/Na);                                                   % number of points in the DFT
noverlap=round(nfft/2);                                             % number of points in the overlap
Sffe=cpsd(fes,fes,hann(nfft),noverlap,nfft,fs);                     % PSD of force applied to the structure
Sff=cpsd(fe,fe,hann(nfft),noverlap,nfft,fs);                        % PSD of force applied to the structure

for n=1:14;
Sww(n,:)=cpsd(dis(n,:),dis(n,:),hann(nfft),noverlap,nfft,fs);       % PSD of displacement response
He(n,:)=tfestimate(fes,dis(n,:),hann(nfft),noverlap,nfft,fs);       % FRF
CoHe(n,:)=mscohere(fes,dis(n,:),hann(nfft),noverlap,nfft,fs);       % coherence
end

dff=1/(nfft*dt);                                                    % frequency resolution
ff = 0:dff:fs/2;                                                    % frequency vector
end


