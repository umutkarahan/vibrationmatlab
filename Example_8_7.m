clear all

%% parameters
E=69e9;                                       % Youngs modulus of aluminium (N/m^2)
rho=2700;                                     % density (kg/m^3);
l=1;b=0.02;d=0.01;S=b*d;I=b*d^3/12;          % geometrical parameters 
z=0.01;n=2*z;                                 % damping ratio and loss factor
Ed=E*(1+j*n);                                 % complex Young's modulus
m=rho*S*l;                                    % mass of the rod% Youngs modulus of aluminium (N/m^2)

%% modal solution
fs=2000;df=0.001;dt=1/fs;                      % frequency parameters
f=0.001:0.01:fs/2;                             % frequency vector
w=2*pi*f;
k=w.^0.5*(rho*S/(Ed*I))^.25;                   % wavenumber

nmax=3;                                       % number of modes
kl(1)=1.87510;kl(2)=4.69409;kl(3)=7.85476;     % kl values 1-3
kl(4)=10.9956;kl(5)=14.1372;                   % kl values 4,5
n=6:nmax;                                      
kl(n)=(2*n-1)*pi/2;                            % kl values > 5

for n=1:nmax
A=(sinh(kl(n))-sin(kl(n)))./(cosh(kl(n))+cos(kl(n)));
x=0.2;                                         % force position
phi1=cosh(kl(n)*x/l)-cos(kl(n)*x/l)-...
A.*(sinh(kl(n)*x/l)-sin(kl(n)*x/l));
x=l;                                           % displacement position 
phi2=cosh(kl(n)*x/l)-cos(kl(n)*x/l)-...
A.*(sinh(kl(n)*x/l)-sin(kl(n)*x/l));
wn=sqrt((E*I)./(rho*S))*(kl(n)).^2;            % natural frequency 
 
Ht(n,:)=phi1*phi2./(m*(wn^2-w.^2+j*2*w*wn*z)); % FRF of each mode
end
Htt=sum(Ht);                                   % overall FRF

%% IRFs
Htd=[Htt fliplr(conj(Htt))];                   % form the double-sided spectrum  
Hm=Htd(1:length(Htd)-1);                       % set the length of the FRF
h=fs*ifft(Hm);                                 % calculation of the IRF
h=circshift(h,10);                             % shift the end of the IRF to the beginning
t=0:dt:(length(h)-1)*dt;                       % time vector

%% plot the results
figure                                              %FRF
semilogx(f,20*log10(abs(Htt)),'linewidth',3,'color',[0.5 0.5 0.5])
set(gca,'fontsize',16)
axis square; grid; axis([1,1010,-150,-30])
xlabel('frequency (Hz)')
ylabel('|FRF| (dB ref 1m/N)')

figure                                              % IRF
plot(t,h,'linewidth',3,'color',[0.5 0.5 0.5])
set(gca,'fontsize',16)
axis square; grid; axis([0,1,-0.02,0.02])
xlabel('time (s)')
ylabel('IRF (m/Ns)')