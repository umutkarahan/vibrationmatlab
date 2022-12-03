%% Combining the measurements of the host structure and the asborber and checking with theory
Example_9_5a                               % delete the figures
Example_9_5b                               % delete the figures and remove clear all at the beginning of the program
%%
Ha=(He'./(-(2*pi*ff).^2))';                 % receptance of SDOF host structure without absorber
Htae=(1./(Kae'+1./He')./(-(2*pi*ff).^2))';  % receptance of SDOF host structure with absorber
%% check on the sum of the measured quantities

mu=0.05;                                    % mass ratio
ma=mu*m;                                    % mass of absorber
wa=wn/(1+mu);                               % tuned frequency of the absorber
ka=wa^2*ma;                                 % stiffness of the absorber
za=sqrt(3/8*mu/(1+mu)^3);                   % optimum absorber damping ratio
ca=2*za*sqrt(ma*ka);                        % optimum absorber damping coefficient
n=0;df=0.01;                                % frequency resolution     
for fc=0:df:100                             % excitation frequency                           
 n=n+1;
 wc=2*pi*fc;
 M=[m 0; 0 ma]; K=[k+ka -ka; -ka ka]; C=[c+ca -ca; -ca ca]; % mass, stiffness and damping matrices
 F=[1;0];                                   % excitation force vector
 D=K-wc.^2*M+j*wc*C;                        % dynamic stiffness of complete system
 HH=inv(D)*F;                               % receptance vector
 Hc(n)=HH(1);                               % receptance of host-structure
end
fc=0:df:100;                                % frequency vector

%% IRF
% calculate IRF of host structure alone
Ha(1)=Ha(2);                                % remove infinity at zero Hz
Hd=[Ha' fliplr(conj(Ha'))]';                % form the double-sided spectrum  
Ht=Hd(1:length(Hd)-1);                      % set the length of the FRF
ha=fs*ifft(Ht);                       % calculation of the IRF

% calculate IRF of host structure with absorber
Htae(1)=Htae(2);                            % remove infinity at zero Hz
Hde=[Htae' fliplr(conj(Htae'))]';           % form the double-sided spectrum  
Hta=Hde(1:length(Hde)-1);                   % set the length of the FRF
hta=fs*ifft(Hta);                     % calculation of the IRF

TT=1/ff(2);tt=0:dt:TT;                      % time vector

%%
figure                                              % FRF modulus
semilogx(ff,20*log10(abs(Ha)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,20*log10(abs(H)),':k','linewidth',4)
hold on
semilogx(ff,20*log10(abs(Htae)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(fc,20*log10(abs(Hc)),'--k','linewidth',4)
axis square,axis([1,100,-120,-40]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');

figure                                                     % phase
semilogx(ff,180/pi*unwrap(angle(Ha)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(f,180/pi*unwrap(angle(H)),':k','linewidth',4)
hold on
semilogx(ff,180/pi*unwrap(angle(Htae)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(fc,180/pi*unwrap(angle(Hc)),'--k','linewidth',4)
axis square,axis([1,100,-200,0]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                           % IRF
plot(tt,ha,'linewidth',2,'color',[0.6 0.6 0.6])
hold on
plot(tt,hta,'k','linewidth',2)
axis square,axis([0,3,-12e-3,12e-3]),grid
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('IRF (m/Ns)');