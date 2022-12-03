clear all 

%% Parameters
m = 1;                              % see Matlab Example 3.1           
k = 500000;%k=10000;                % large stiffness to illustrate aliasing at low frequencies
z = 0.01; c = 2*z*sqrt(m*k);        %  and a small stiffness to illustrate aliasing at high frequencies
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;   

%% Time and frequency parameters
T=10;                               % see Matlab Example 3.1
fs=400;                             
dt=1/fs; t=0:dt:T;                  
df=1/T; f=0:df:fs;                  
N=length(t);                              

%% Theoretical receptance IRF
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % IRF

%% Calculation of DFT
H=dt*fft(h);                        % Calculation of receptance FRF

%% Theoretical FRF
dff=df; fr=0:dff:fs;                % frequency resolution/frequency vector
w=2*pi*fr;                          % angular frequency
HH=1./(k-w.^2*m+j*w*c);             % theoretical FRF

%% Calculation of aliased response
for p=1:20;
f1=(p-1)*fs:dff:p*fs;               % frequency vector        
w1=2*pi*f1;
HP(p,:)=1./(k-w1.^2*m+j*w1*c);      % aliased FRF for +ve frequency
f2=-p*fs:dff:-(p-1)*fs;             % frequency vector 
w2=2*pi*f2;
HM(p,:)=1./(k-w2.^2*m+j*w2*c);      % aliased FRF for -ve frequency
end
MP=sum(HP); MS=(sum(HM));           % summing all aliases
HT1=MP+MS;                          % total aliased FRF
HT=HT1(1:(N+1)/2);
HA=[HT fliplr(conj(HT))];           % double sided spectrum
HA1=HA(1:length(f));

H0=1/k-1/(12*fs^2*m);               % value of aliased FRF at f=0
H1_2=1/(4*fs^2*m);                  % value of aliased FRF at f=fs/2

%% Plot results
figure (1)
semilogx(fr,20*log10(abs(HH)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,20*log10(abs(H)),'Color',[.6 .6 .6],'linewidth',4)
hold on
semilogx(fr,20*log10(HT1),'--k','linewidth',3)
hold on
% plot(fs/2,20*log10(abs(1/((pi*fs).^2*m))),'--ok','linewidth',2,'Markersize',8)
% plot(fs/2,20*log10(H1_2),'--ok','linewidth',2,'Markersize',8)
axis([1,fs/2,-140,-40])
plot(1,20*log10(H0),'ok','linewidth',2,'Markersize',8)
plot(1,20*log10(abs(1/k)),'ok','linewidth',2,'Markersize',8)
axis([1,fs/2,-130,-80])
xticks([1e0 1e1 1e2])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
grid;axis square

%%
figure (2)
semilogx(fr,180/pi*(angle(HH)),'linewidth',4,'Color',[.7 .7 .7])
hold on
semilogx(f,180/pi*unwrap(angle(H)),'--k','linewidth',2)
axis square
axis([1,fs/2,-200,0])
set(gca,'fontname', 'arial','fontsize',24)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
