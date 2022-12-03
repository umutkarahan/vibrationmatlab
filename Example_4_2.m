clear all 

%% Parameters
m = 1;                              % see Matlab Example 3.1             
k = 10000;                            
z = 0.01; c = 2*z*sqrt(m*k);         
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    

%% Time and frequency parameters
T=10;                              % see Matlab Example 3.1
fs=1000;                            
dt=1/fs; t=0:dt:T;                  
df=1/T; f=0:df:fs;                  
N=length(t);                              

%% Theoretical velocity IRF
hv1= 1/(m*wd)*exp(-z*wn*t);
hv2= (wd*cos(wd*t)-z*wn*sin(wd*t));
hv= hv1.*hv2;                       % IRF

%% Calculation of DFT
Hv=dt*fft(hv);                      % calculation of mobility FRF

%% Theoretical mobility FRF
dff=df; fr=0:dff:fs;                % see Matlab example 3.1
w=2*pi*fr;                         
HV=j*w./(k-w.^2*m+j*w*c);           % mobility FRF

%% Calculation of aliased response
for p=1:20;
 f1=(p-1)*fs:dff:p*fs;                                        % frequency vector
 w1=2*pi*f1;
 HP(p,:)=j*w1./(k-w1.^2*m+j*w1*c);                            % aliased FRF for +ve frequencies
 EP(p,:)=j./(w1*m).*(1-exp(j*w1*dt/2).*sinc(w1/(2*pi)*dt));   % aliased additional component for +ve frequencies
 f2=-p*fs:dff:-(p-1)*fs;                                      % frequency vector
 w2=2*pi*f2;
 HM(p,:)=j*w2./(k-w2.^2*m+j*w2*c);                            % aliased FRF for -ve frequencies
 EM(p,:)=j./(w2*m).*(1-exp(j*w2*dt/2).*sinc(w2/(2*pi)*dt));   % aliased additional component for -ve frequencies
end
MP=sum(HP+EP); MS=(sum(HM+EM));                               % sum of all aliases
HT1=MP+MS;                                                    % total aliased FRF
HT=HT1(1:(N+1)/2);
HA=[HT fliplr(conj(HT))];                                     % double sided spectrum
HA1=HA(1:length(f));

H0=1/(2*fs*m);                                                % value of aliased FRF at f=0
H1_2=H0;                                                      % value of aliased FRF at f=fs/2
%% Plot the results
figure (1)                                                    % modulus
semilogx(fr,20*log10(abs(HV)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,20*log10(abs(Hv)),'Color',[.6 .6 .6],'linewidth',4)
hold on
semilogx(fr,20*log10(abs(HA1)),'--k','linewidth',3)
hold on
plot(fs/2,20*log10(abs(1/(pi*fs*m))),'ok','linewidth',2,'Markersize',8)
plot(fs/2,20*log10(H1_2),'ok','linewidth',2,'Markersize',8)
plot(0.1,20*log10(H0),'ok','linewidth',2,'Markersize',8)
axis([0.1,fs/2,-90,0])
xticks([1e0 1e2])
set(gca,'fontsize',16)
%%
figure (2)                                                    % phase
xlabel('frequency (Hz)');
ylabel('|mobility| (dB ref 1 m/Ns)');
grid;axis square
semilogx(fr,180/pi*(angle(HV)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(fr,180/pi*(angle(HT1)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,180/pi*unwrap(angle(Hv)),'--k','linewidth',3)
axis square
axis([0.1,fs/2,-90,90])
xticks([1e0 1e2])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
grid;axis square