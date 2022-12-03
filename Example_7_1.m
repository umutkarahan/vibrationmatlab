%% calculation of the FRF due to a force impulse and a chirp
clear all 
%% Parameters
m = 1;                              % mass           
k=10000;                            % stifness  
z = 0.1; c = 2*z*sqrt(m*k);         % damping 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % natural frequency
fn=wn/(2*pi);                       % natural frequency
Tn=1/fn;                            % natural period 
SNRf=20;SNRx=100;                   % added noise SNRs
%% half-sine pulse
Tc=Tn*0.02;                              % contact time for half sine pulse
fs=5000;                                 % Sampling frequency
dt=1/fs; t1=0:dt:Tc;                     % time vector
fp = sin(pi*t1/Tc);                      % half sine pulse
Nz=600*length(t1);                       % Number of zeros to add
f = [zeros(Nz/10,1)',fp zeros(Nz,1)'];   % zero-padded force signal
t = 0:dt:(length(f)-1)*dt;               % time vector for extended signal
N=length(t);TT=N*dt;Tm=max(t);
[xc]=calculate(f,m,k,c,wd,wn,z,t,fs);    % calculate displacement
[fwn,xwn,Sff,Sxx,H1,H1a,H2,H2a,coh]=FRF(f,xc,dt,Tm,SNRf,SNRx); % calculate frequency domain quantities
Sffs=Sff;Sxxs=Sxx;H1s=H1;H1as=H1a;H2s=H2;H2as=H2a;cohs=coh;    % frequency domain quantities
figure
plots(t,fwn,TT,xwn)        
%% chirp
T=TT/2;tt=0:dt:T;                        % time vector  
f1=1;f2=200;                             % upper and lower frequencies
a=2*pi*(f2-f1)/(2*T); b=2*pi*f1;         % coefficients
fc=sin(a*tt.^2+b*tt);                    % chirp signal
f=[fc zeros(1,length(f)-length(fc))];    % zero padded force signal
[xc]=calculate(f,m,k,c,wd,wn,z,t,fs);    % calculate displacement
[fwn,xwn,Sff,Sxx,H1,H1a,H2,H2a,coh]=FRF(f,xc,dt,Tm,SNRf,SNRx); % calculate frequency domain properties
Sffc=Sff;Sxxc=Sxx;H1c=H1;H1ac=H1a;H2c=H2;H2ac=H2a;cohc=coh;
figure
plots(t,fwn,TT,xwn)

%% theoretical FRF
df=1/(N*dt);ff=0:df:fs-df;               % frequency vector
w=2*pi*ff;
Ht=1./(k-w.^2*m+j*w*c);                  % theoretical FRF

%% Frequency domain plots
figure                                                 % plot of force PSD
semilogx(ff,10*log10(mean(Sffc)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,10*log10(mean(Sffs)),'linewidth',4,'color',[0.2 0.2 0.2])
axis square; grid; axis([1,1000,-80,-30])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('force PSD (dB ref 1 N^2/Hz)');

figure                                                 % plot of displacement PSD
semilogx(ff,10*log10(mean(Sxxc)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,10*log10(mean(Sxxs)),'linewidth',4,'color',[0.2 0.2 0.2])
axis square; grid; axis([1,1000,-200,-80])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('displ. PSD (dB ref 1 m^2/Hz)');

figure                                                 % plot of H1 estimator
semilogx(ff,20*log10(abs(H1a)),'linewidth',4,'color',[0.7 0.7 0.7])
hold on
semilogx(ff,20*log10(abs(H1c)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,20*log10(abs(H1s)),'linewidth',4,'color',[0.2 0.2 0.2])
hold on
semilogx(ff,20*log10(abs(Ht)),':k','linewidth',4)
axis square; grid; axis([1,1000,-120,-60])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('displ./force (dB ref 1 m/N)');

figure                                                % plot of H2 estimator
semilogx(ff,20*log10(abs(H1ac)),'linewidth',4,'color',[0.7 0.7 0.7])
hold on
semilogx(ff,20*log10(abs(H2c)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,20*log10(abs(H2s)),'linewidth',4,'color',[0.2 0.2 0.2])
hold on
semilogx(ff,20*log10(abs(Ht)),':k','linewidth',4)
axis square; grid; axis([1,1000,-120,-60])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('displ./force (dB ref 1 m/N)');

figure                                              % plot of phase
semilogx(ff,180/pi*angle(H1a),'linewidth',4,'color',[0.7 0.7 0.7])
hold on
semilogx(ff,180/pi*unwrap(angle(H1c)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,180/pi*unwrap(angle(H1s)),'linewidth',4,'color',[0.2 0.2 0.2])
hold on
semilogx(ff,180/pi*unwrap(angle(Ht)),':k','linewidth',4)
axis square; grid; axis([1,1000,-200,200])
xticks([1 10 100 1000])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                              % plot of coherence
semilogx(ff,cohc,'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(ff,cohs,'linewidth',4,'color',[0.2 0.2 0.2])
axis square; grid; axis([1,1000,0,1])
xticks([1 10 100 1000])
yticks([0 0.2 0.4 0.6 0.8 1])
set(gca,'fontsize',20)
xlabel('frequency (Hz)');
ylabel('coherence');

%% Functions
function [xc]=calculate(f,m,k,c,wd,wn,z,t,fs)      % calculate displacement
%% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % IRF
%% Convolution
xc=conv(h,f)/fs;                                   % displacement response   
xc=xc(1:length(f));                 
end

function [fwn,xwn,Sff,Sxx,H1,H1a,H2,H2a,coh]=FRF(f,xc,dt,Tm,SNRf,SNRx) % calculate frequency domain properties
for n=1:16
fwn=awgn(f,SNRf,'measured','dB');                  % add random noise
xwn=awgn(xc,SNRx,'measured','dB');                 % add random noise

F=fft(fwn)*dt;                                     % fft of force
X=fft(xwn)*dt;                                     % fft of displacement
Sff(n,:)=F.*conj(F)/Tm;                            % force PSD                     
Sxx(n,:)=X.*conj(X)/Tm;                            % displacement PSD             
Sfx(n,:)=X.*conj(F)/Tm;                            % CPSD fx
Sxf(n,:)=F.*conj(X)/Tm;                            % CPSD xf
end

coh=abs(mean(Sxf)).^2./(mean(Sxx).*mean(Sff));     % coherence
H1=mean(Sfx)./mean(Sff); H1a=Sfx./Sff;             % H1 estimator
H2=mean(Sxx)./mean(Sxf); H2a=Sxx./Sxf;             % H2 estimator
end

function plots(t,fwn,TT,xwn)
subplot(2,2,1)
plot(t,fwn,'linewidth',3,'Color',[.6 .6 .6]),grid
set(gca,'fontsize',16)
axis([0,TT,1.1*min(fwn),1.1*max(fwn)])
xlabel('time (s)');
ylabel('force (N)');
subplot(2,2,3)
plot(t,xwn,'linewidth',3,'Color',[.6 .6 .6]),grid
set(gca,'fontsize',16)
axis([0,TT,1.1*min(xwn),1.1*max(xwn)])
xlabel('time (s)');
ylabel('displacement (m)');
 end
