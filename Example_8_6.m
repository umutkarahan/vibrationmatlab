clear all

%% parameters
E=69e9;                                       % Youngs modulus of aluminium (N/m^2)
rho=2700;                                     % density (kg/m^3);
l=10;b=0.02;d=0.01;S=b*d;                     % geometrical parameters 
z=0.01;n=2*z;                                 % damping ratio and loss factor
Ed=E*(1+j*n);                                 % complex Young's modulus
m=rho*S*l;                                    % mass of the rod

%% modal solution
fs=20000;df=0.001;dt=1/fs;                     % frequency parameters
f=0.001:0.01:fs/2;                             % frequency vector
w=2*pi*f;
for n=1:40                                         % 40 modes (can change this number)
    x=l;                                           % force position                                      
    phi1=sqrt(2)*sin((n-1/2)*pi*x/l);              % mode shape at force position
    x=0.5*l;                                       % response position 
    phi2=sqrt(2)*sin((n-1/2)*pi*x/l);              % mode shape at response position
    wn=(n-1/2)*pi/l*sqrt(E/rho);                   % natural frequencies
    Ht(n,:)=phi1*phi2./(m*(wn^2-w.^2+j*2*w*wn*z)); % FRF for each mode
end
Htt=sum(Ht);                                   % overall FRF

%% IRF
Htd=[Htt fliplr(conj(Htt))];                   % form the double-sided spectrum  
Hm=Htd(1:length(Htd)-1);                       % set the length of the FRF
h=fs*ifft(Hm);                                 % calculation of the IRF
h=circshift(h,100);                            % shift the end of the IRF to the beginning

t=0:dt:(length(h)-1)*dt;                       % time vector

%% plot the results
figure                                              % FRF
semilogx(f,20*log10(abs(Htt)),'linewidth',3,'color',[0.5 0.5 0.5])
set(gca,'fontsize',16)
axis square; grid; axis([10,10000,-180,-90])
xlabel('frequency (Hz)')
ylabel('|FRF| (dB ref 1m/N)')

figure                                              % IRF
plot(t,h,'linewidth',3,'color',[0.5 0.5 0.5])
set(gca,'fontsize',16)
axis square; grid; axis([0,0.1,-4e-4,4e-4])
xlabel('time (s)')
ylabel('IRF (m/Ns)')