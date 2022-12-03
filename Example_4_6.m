clear all

%% Parameters
m = 1;                              % see Matlab example 3.1           
k = 10000;                            
z = 0.1; c = 2*z*sqrt(m*k);          
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    

%% Time and frequency parameters
fs=200;                             % see Matlab example 3.1
T=1;                                
dt=1/fs; t=0:dt:T;                  
df=1/T; f=0:df:fs;                  

%% Theoretical IRFs
%displacement
hdt= 1/(m*wd)*exp(-z*wn*t).*sin(wd*t);         % displacement IRF
%velocity
hv1= 1/(m*wd)*exp(-z*wn*t);
hv2= (wd*cos(wd*t)-z*wn*sin(wd*t));
hvt= hv1.*hv2;                                 % velocity IRF
%acceleration
ha1= 1/(m*wd)*exp(-z*wn*t);
ha2=-(wd^2*sin(wd*t) + z*wn*wd*cos(wd*t));
ha3= -z*wn*(wd*cos(wd*t) - z*wn*sin(wd*t));
hat=ha1.*(ha2+ha3);hat(1)=1/(dt*m)+hat(1);     % acceleration IRF

%% Calculation of DFTs
Hd=dt*fft(hdt);                                % receptance
Hv=dt*fft(hvt);                                % mobility
Ha=dt*fft(hat);                                % accelerance

%% Theoretical FRFs
fr=0:df:fs/2;                                  % frequency vector                 
w=2*pi*fr;                                     
HD=1./(k-w.^2*m+j*w*c);                        % receptance
HV=j*w.*HD;                                    % mobility
HA=j*w.*HV;                                    % accelerance

%% Form double-sided spectra
Hda=[HD fliplr(conj(HD))];HDd=Hda(1:length(t)); % receptance 
Hva=[HV fliplr(conj(HV))];HVd=Hva(1:length(t)); % mobility
Haa=[HA fliplr(conj(HA))];HAd=Haa(1:length(t)); % accelerance   
%% inverse Fourier transforms
hdi=fs*ifft(HDd);                              % displacement IRF
hvi=fs*ifft(HVd);                              % velocity IRF
hai=fs*ifft(HAd);                              % acceleration IRF
%% calculation of convolved sinc function with IRF
a=0.05;dtt=a*dt;tt=0:dtt:T;
ts=-T/2:dtt:T/2;
W=fs*sinc(fs*ts);
%displacement
hd= 1/(m*wd)*exp(-z*wn*tt).*sin(wd*tt);        % displacement IRF
hdc=conv(hd,W)*dtt;hdcc=hdc(1:length(tt));     % convolve IRF with sinc function
%velocity
hv1= 1/(m*wd)*exp(-z*wn*tt);                   
hv2= (wd*cos(wd*tt)-z*wn*sin(wd*tt));
hv=hv1.*hv2;                                   % velocity IRF 
hvc=conv(hv,W)*dtt;hvcc=hvc(1:length(tt));     % convolve IRF with sinc function
%acceleration
ha1= 1/(m*wd)*exp(-z*wn*tt);
ha2=-(wd^2*sin(wd*tt) + z*wn*wd*cos(wd*tt));
ha3= -z*wn*(wd*cos(wd*tt) - z*wn*sin(wd*tt));
ha=ha1.*(ha2+ha3); ha(1)=1/(dtt*m)+ha(1);      % acceleration IRF
hac=conv(ha,W)*dtt;hacc=hac(1:length(tt));     % convolve IRF with sinc function
%% shifting the IRFs to seee the acausality
hdcirc=circshift(hdt',100);
hdicirc=circshift(hdi',100);
hdcccirc=circshift(hdcc',(length(tt)+1)/2+100/a);
hvcirc=circshift(hvt',100);
hvicirc=circshift(hvi',100);
hvcccirc=circshift(hvcc',(length(tt)+1)/2+100/a);
hacirc=circshift(hat',100);
haicirc=circshift(hai',100);
hacccirc=circshift(hacc',(length(tt)+1)/2+100/a);
%% plot the results
%displacement
figure (1)                                                  % modulus
plot(f,20*log10(abs(Hd)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,20*log10(abs(HDd)),'ok','Linewidth',2,'MarkerSize',8)
axis([0,fs,-120,-60])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
grid;axis square
%
figure (2)                                                   % phase
plot(f,180/pi*(angle(Hd)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,180/pi*(angle(HDd)),'ok','Linewidth',2,'MarkerSize',8)
axis([0.1,fs,-200,200])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
grid;axis square

figure (3)                                                   % IRF
plot(tt,hdcccirc,'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(t,hdicirc,'ok','Linewidth',2,'MarkerSize',8)
hold on
plot(t,hdcirc,'k','Markersize',6)
axis([0.4,1,-0.01,0.01])
grid;axis square
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time(s)');
ylabel('displacement IRF (m/Ns)');

%velocity
figure (4)                                                   % modulus
plot(f,20*log10(abs(Hv)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,20*log10(abs(HVd)),'ok','Linewidth',2,'MarkerSize',8)
axis([0,fs,-70,-20])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('|mobility| (dB ref 1 m/Ns)');
grid;axis square
%
figure (5)                                                   % phase
plot(f,180/pi*(angle(Hv)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,180/pi*(angle(HVd)),'ok','Linewidth',2,'MarkerSize',8)
grid;axis square
axis([0,fs,-100,100])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure (6)                                                   % IRF
plot(tt,hvcccirc,'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(t,hvicirc,'ok','linewidth',2)
hold on
plot(t,hvcirc,'k','Linewidth',2)
axis([0.4,1,-1,1.2])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time(s)');
ylabel('velocity IRF (m/Ns^2)');
grid;axis square

figure (7)                                                   % modulus
plot(f,20*log10(abs(Ha)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,20*log10(abs(HAd)),'ok','Linewidth',2,'MarkerSize',8)
axis([0,fs,-40,20])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('|accelerance| (dB ref 1 m/Ns^2)');
grid;axis square
%
figure (8)                                                   % phase
plot(f,180/pi*unwrap(angle(Ha)),'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(f,180/pi*(angle(HAd)),'ok','Linewidth',2,'MarkerSize',8)
grid;axis square
axis([0,fs,-200,200])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure (9)                                                   % IRF                         
plot(tt,hacccirc,'Linewidth',4,'Color',[.6 .6 .6])
hold on
plot(t,haicirc,'ok','Linewidth',2,'MarkerSize',8)
hold on
plot(t,hacirc,'k','Linewidth',2)
axis([0.4,1,-150,200])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('time(s)');
ylabel('acceleration IRF (m/Ns^3)');
grid;axis square