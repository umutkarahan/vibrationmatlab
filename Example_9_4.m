clear all

Example_9_3                            % run the program for the last example and delete the figures

%% mode shapes
xx=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75];
% obtain the values of the imaginary part of the receptance at the first 3 resonance frequencies
p1=[0.13 0.29 0.49 0.74 1.03 1.35 1.7 2.07 2.45 2.85 3.25 3.67 4.07 4.48];
p2=[0.11 0.22 0.34 0.44 0.50 0.53 0.51 0.44 0.31 0.15 -0.05 -0.27 -0.50 -0.74];
p3=[0.08 0.13 0.16 0.16 0.12 0.04 -0.04 -0.1 -0.14 -0.14 -0.09 0 0.1 0.22];

% calculate the theoretical mode shapes for the first 3 resonance frequencies
kl=1.87510;
[xf,pt]=modeCalc(kl,l);                       % function for mode shapes
p1t=pt;
kl=4.69409;
[xf,pt]=modeCalc(kl,l);                       % function for mode shapes
p2t=pt;
kl=7.85476;
[xf,pt]=modeCalc(kl,l);                       % function for mode shapes
p3t=pt;

%% modal properties
% construct the modal model W1/Fe; 
f1=8.13;f2=51.2;f3=143.3;                                               % estimate the natural frequencies from receptance FRF
x1=0.1315e-3;x2=0.1138e-3;x3=7.83e-5;                                   % estimate the responses at the resonance frequencies from the Imag. part of the receptance
z1=(8.27-8.07)/(2*f1); z2=(51.73-50.73)/(2*f2);z3=(144.8-141.9)/(2*f3); % estimate the modal damping ratios using the half-power points
w=2*pi*ff;w1=2*pi*f1;w2=2*pi*f2;w3=2*pi*f3;                             % frequency and natural frequencies in rad/s
A1=x1.*(2*z1*w1^2);A2=x2.*(2*z2*w2^2);A3=x3.*(2*z3*w3^2);               % determine the modal constants

Rm(1,:)=A1./(w1^2-w.^2+j*2*z1*w*w1)+A2./(w2^2-w.^2+j*2*z2*w*w2)+A3./(w3^2-w.^2+j*2*z3*w*w3); % estimated receptance from modal parameters 
AA=1/(He(1,1)-Rm(1));                                                   % estimate the residual
Rmm(1,:)= Rm(1,:) + w./w*1./AA;                                         % add the residual

% construct the modal model W7/Fe;
f1=8.13;f2=51.2;f3=143.3;                                               % estimate the natural frequencies from receptance FRF
x1=1.699e-3;x2=0.5085e-3;x3=-3.55e-5;                                   % estimate the responses at the resonance frequencies from the Imag. part of the receptance
z1=(8.27-8.07)/(2*f1); z2=(51.73-50.73)/(2*f2);z3=(144.8-141.9)/(2*f3); % estimate the modal damping ratios using the half-power points
w=2*pi*ff;w1=2*pi*f1;w2=2*pi*f2;w3=2*pi*f3;                             % frequency and natural frequencies in rad/s
A1=x1.*(2*z1*w1^2);A2=x2.*(2*z2*w2^2);A3=x3.*(2*z3*w3^2);               % determine the modal constants

Rm(7,:)=A1./(w1^2-w.^2+j*2*z1*w*w1)+A2./(w2^2-w.^2+j*2*z2*w*w2)+A3./(w3^2-w.^2+j*2*z3*w*w3); % estimated receptance from modal parameters 
AA=1/(He(7,1)-Rm(7,1)),                                                 % estimate the residual
Rmm(7,:)= Rm(7,:) + w./w*1./AA;                                         % add the residual

% construct the modal model W14/Fe;
f1=8.13;f2=51.2;f3=143.3;                                               % estimate the natural frequencies from receptance FRF
x1=4.483e-3;x2=-0.7371e-3;x3=0.2183e-3;                                 % estimate the responses at the resonance frequencies from the Imag. part of the receptance
z1=(8.27-8.07)/(2*f1); z2=(51.73-50.73)/(2*f2);z3=(144.8-141.9)/(2*f3); % estimate the modal damping ratios using the half-power points
w=2*pi*ff;w1=2*pi*f1;w2=2*pi*f2;w3=2*pi*f3;                             % frequency and natural frequencies in rad/s
A1=x1.*(2*z1*w1^2);A2=x2.*(2*z2*w2^2);A3=x3.*(2*z3*w3^2);               % determine the modal constants

Rm(14,:)=A1./(w1^2-w.^2+j*2*z1*w*w1)+A2./(w2^2-w.^2+j*2*z2*w*w2)+A3./(w3^2-w.^2+j*2*z3*w*w3); % estimated receptance from modal parameters 
AA=1/(He(14,1)-Rm(14,1));                                               % estimate the residual
Rmm(14,:)= Rm(14,:) + 0*w./w*1./AA;                                       % add the residual
%% plot the results
figure                                     % first mode-shape
plot(xf,p1t/max(p1t),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(xx,p1/4.48,'ok','Markersize',10,'linewidth',2)
axis square,axis([0,0.75,-0.1,1]),grid
set(gca,'fontsize',16)
xlabel('beam position (m)');ylabel('mode shape');

figure                                      % second mode-shape
plot(xf,-p2t/min(p2t),'linewidth',4,'color',[0.6 0.6 0.6])  
hold on
plot(xx,p2/0.74,'ok','Markersize',10,'linewidth',2)
axis square,axis([0,0.75,-1,1]),grid
set(gca,'fontsize',16)
xlabel('beam position (m)');ylabel('mode shape');

figure                                     % third mode-shape
plot(xf,p3t/max(p3t),'linewidth',4,'color',[0.6 0.6 0.6])              
hold on
plot(xx,p3/0.22,'ok','Markersize',10,'linewidth',2)
axis square,axis([0,0.75,-1,1]),grid
set(gca,'fontsize',16)
xlabel('beam position (m)');ylabel('mode shape');
%%
p=14;                                                    % response position
figure                                                   % modulus of FRF
semilogx(ff,20*log10(abs(He(p,:))),'linewidth',4,'color',[0.6 0.6 0.6])
hold on 
semilogx(ff,20*log10(abs(Rm(p,:))),':k','linewidth',4)
hold on 
semilogx(ff,20*log10(abs(Rmm(p,:))),'--k','linewidth',4)
axis square,axis([1,1000,-160,-40]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');

figure                                                    % phase
semilogx(ff,180/pi*(angle(He(p,:))),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
semilogx(ff,180/pi*(angle(Rm(p,:))),':k','linewidth',4)
hold on
semilogx(ff,180/pi*(angle(Rmm(p,:))),'--k','linewidth',4)
axis square,axis([1,1000,-200,+200]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

%%
function [xf,pt]=modeCalc(kl,l);                                         % function for mode shape calculation
A=(sinh(kl)-sin(kl))./(cosh(kl)+cos(kl));
xf=0:0.01:l;                                         
pt=cosh(kl*xf/l)-cos(kl*xf/l)-...
A.*(sinh(kl*xf/l)-sin(kl*xf/l));
end



