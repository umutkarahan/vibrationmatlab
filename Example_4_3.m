clear all 

%% Parameters
m = 10;                              % see Matlab Example 3.1           
k = 10000;                             
z = 0.01; c = 2*z*sqrt(m*k);         % choose z=0.001 for light damping
wn = sqrt(k/m); wd = sqrt(1-z^2)*wn;     
fn = sqrt(k/m)/(2*pi);                 

%% Time and frequency parameters
T=100;T=500;                        % see Matlab example 3.1
fs=1000;                            % choose T=500 for low damping to avoid data Truncation
dt=1/fs; t=0:dt:T;                  
df=1/T; f=0:df:fs;                  
N=length(t);                              

%% Theoretical acceleration IRF
ha1= 1/(m*wd)*exp(-z*wn*t);
ha2=-(wd^2*sin(wd*t) + z*wn*wd*cos(wd*t));
ha3= -z*wn*(wd*cos(wd*t) - z*wn*sin(wd*t));
ha=ha1.*(ha2+ha3);
ha(1)=1/(dt*m)+ha(1);               % IRF including the delta function and oscillatory term


%% Calculation of DFT
Ha=dt*fft(ha);                      % calculation of accelerance FRF

%% Theoretical FRF
dff=df; fr=0:dff:fs;                % see Matlab example 3.1
w=2*pi*fr;                          
HA=-w.^2./(k-w.^2*m+j*w*c);         

%% Calculation of aliased response
for p=1:20;
 f1=(p-1)*fs:dff:p*fs;                            % frequency vector
 w1=2*pi*f1;
 HP(p,:)=-w1.^2./(k-w1.^2*m+j*w1*c)-1/m;          % aliased FRF for +ve frequencies
 AP(p,:)=-j*2*z*wn./(m*w1).*(1-exp(j*w1*dt/2)).*sinc(w1/(2*pi)*dt); 
 DP(p,:)=1/m*(sinc(w1/(2*pi)*dt)).^2;             % aliased additional component                  
 f2=-p*fs:dff:-(p-1)*fs;
 w2=2*pi*f2;
 HM(p,:)=-w2.^2./(k-w2.^2*m+j*w2*c)-1/m;          % aliased FRF for -ve frequencies
 AM(p,:)=-j*2*z*wn./(m*w2).*(1-exp(j*w2*dt/2)).*sinc(w2/(2*pi)*dt);
 DM(p,:)=1/m*(sinc(w2/(2*pi)*dt)).^2;             % aliased additional component
end
MP=sum(HP+AP+DP);MS=sum(HM+AM+DM);                % sum all aliases
HT1=MP+MS;                                        % total aliased FRF
HT=HT1(1:(N+1)/2);
HAA=[HT fliplr(conj(HT))];                        % double side spectrum
HA1=HAA(1:length(f));
H0=(1/m*(pi^2/3*fn^2/fs^2 - 2*pi*z*fn/fs));       % value of aliased FRF at f=0

%% Plot the results
figure (1)                                        % modulus
semilogx(fr,20*log10(abs(HA)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,20*log10(abs(Ha)),'Color',[.6 .6 .6],'linewidth',3)
hold on
semilogx(f,20*log10(abs(HA1)),'--k','linewidth',3)
plot(df,20*log10(abs((2*pi*df)^2/k)),'ok','linewidth',2,'Markersize',8)
plot(df,20*log10(abs(H0)),'ok','linewidth',2,'Markersize',8)
xticks([1e-2 1e0 1e2])
axis([df,fs/2,-140,20])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|accelerance| (dB ref 1 m/Ns^2)');
grid;axis square

figure (2)                                        % phase  
semilogx(fr,180/pi*(angle(HA)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(fr,360+180/pi*unwrap(angle(HT1)),'linewidth',4,'Color',[.6 .6 .6])
hold on
semilogx(f,180/pi*unwrap(angle(Ha)),'--k','linewidth',3)
xticks([1e-2 1e0 1e2])
grid,axis square
axis([0.01,fs/2,0,200])
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');
