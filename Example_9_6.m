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
f=0:df:fs/2;                                  % frequency vector
w=2*pi*f;

%% positions along the beam
xp=0.1:0.05:l;                                % positions along the beam
fp=1;                                         % force position in displ. vector
ap=14;                                        % absorber position 
op=1;                                         % measurement position 

%% beam FRFs
n=0;
for x=0.1:0.05:l;
n=n+1;
xf=xp(fp);                                        % position of excitation force
[Htt,wnn] = calcFRF(E,I,rho,S,z,m,l,x,xf,w,fs);   % calculate beam FRFs and IRFs
H1(n,:)=Htt;                                      % beam FRFs wrt excitation force
xf=xp(ap);                                        % position of absorber
[Htt,wnn] = calcFRF(E,I,rho,S,z,m,l,x,xf,w,fs);   % calculate beam FRFs and IRFs
H2(n,:)=Htt;                                      % beam FRFs wrt absorber position
end

%% modal properties (first 3 modes)
x=0.1:0.05:l;                                     % position along the beam 
fn=wnn/(2*pi);                                    % natural freqs in Hz
z1=(8.27-8.07)/(2*fn(1)); z2=(51.73-50.73)/(2*fn(2));z3=(144.8-141.9)/(2*fn(3));   % modal damping ratios
xx1=max(abs(H2(ap,700:1000)));                    % max value of 1st modal response
xx2=max(abs(H2(ap,5000:5500)));                   % max value of 2nd modal response
xx3=max(abs(H2(ap,14000:15000)));                 % max value of 3rd modal response
A1=xx1.*(2*z1*wnn(1)^2);                          % 1st modal constant           
A2=xx2.*(2*z2*wnn(2)^2);                          % 2nd modal constant  
A3=xx3.*(2*z3*wnn(3)^2);                          % 3rd modal constant

%% vibration absorber design
mu=0.05;                                          % mass ratio
ma=0.05/A1;                                       % mass of absorber for 1st mode
wa=wnn(1)/(1+mu);                                 % absorber natural frequency for 1st mode
ka=wa^2*ma;                                       % stiffness of the absorber
za=sqrt(3/8*mu/(1+mu)^3);                         % absorber damping ratio
ca=2*za*sqrt(ma*ka);                              % absorber damping coefficient
Ka=-w.^2.*ma.*(ka+j*w*ca)./(ka-w.^2*ma+j*w*ca);   % absorber dynamic stiffness 
Ha=1./Ka;                                         % absorber receptance

%% FRF of beam with absorber attached
H1c=H1(op,:)-H1(ap,:).*H2(op,:)./(H2(ap,:)+Ha);   % Beam FRF with absorber attached

%% IRF
% calculate IRF of host structure alone
Hd=[H1(op,:) fliplr(conj((H1(op,:))))];           % form the double-sided spectrum  
Ht=Hd(1:length(Hd)-1);                            % set the length of the FRF
ha=fs*ifft(Ht);                                   % calculation of the IRF

% calculate IRF of host structure with absorber
Hde=[H1c fliplr(conj(H1c))];                      % form the double-sided spectrum  
Hta=Hde(1:length(Hde)-1);                         % set the length of the FRF
hta=fs*ifft(Hta);                                 % calculation of the IRF

TT=1/f(2);tt=0:dt:TT;                             % time vector

%% plot the results
figure                                                     % FRF modulus
semilogx(f,20*log10(abs(H1(op,:))),'linewidth',4,'color',[0.6 0.6 0.6]);
hold on 
semilogx(f,20*log10(abs(H1c)),'k','linewidth',4);
axis square,axis([1,1000,-160,-40]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');
%% 
figure                                                    % IRF
plot(tt,ha,'linewidth',3,'color',[0.6 0.6 0.6])
hold on
plot(tt,hta,'linewidth',3,'color',[0.3 0.3 0.3])
axis square,axis([0,3,-15e-3,15e-3]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('IRF (m/Ns)');
%%
function [Htt,wnn]=calcFRF(E,I,rho,S,z,m,l,x,xf,w,fs)   % function to calculate FRF and IRF 
nmax=10;                                                % number of modes
kl(1)=1.87510;kl(2)=4.69409;kl(3)=7.85476;              % kl values 1-3
kl(4)=10.9956;kl(5)=14.1372;                            % kl values 4,5
n=6:nmax;                                      
kl(n)=(2*n-1)*pi/2;                                     % kl values > 5
for n=1:nmax;
A=(sinh(kl(n))-sin(kl(n)))./(cosh(kl(n))+cos(kl(n)));                                               
phi1=cosh(kl(n)*xf/l)-cos(kl(n)*xf/l)-...
A.*(sinh(kl(n)*xf/l)-sin(kl(n)*xf/l));

phi2=cosh(kl(n)*x/l)-cos(kl(n)*x/l)-...                 % response position
A.*(sinh(kl(n)*x/l)-sin(kl(n)*x/l));
wn=sqrt((E*I)./(rho*S))*(kl(n)).^2;                     % natural frequency 
wnn(n)=wn;
Ht(n,:)=phi1*phi2./(m*(wn^2-w.^2+j*2*w*wn*z));          % FRF of each mode
end
Htt=sum(Ht);                                            % overall receptance FRF
 end


