%% Development of a modal model for the system in Fig 9.2
clear all
Exercise_9_1                        % delete the figures
%% calculation of receptance
Re=He./(-(2*pi*ff).^2)';            % receptance

%% modal properties
f1=9.84; f2=25.71;                                          % estimation of natural frequencies from receptance FRF
x1=0.0148;x2=0.0003222;                                     % estimation of amplitudes at the resonance frequencies from the Imag. part of the receptance
z1=(9.9-9.78)/(2*9.84);z2=(26.1-25.23)/(2*25.71)            % estimation of the modal damping ratios using the half-power points
w=2*pi*ff;w1=2*pi*f1;w2=2*pi*f2;                            % frequency and natrural frequencies in rad/s
A1=x1.*(2*z1*w1^2); A2=x2.*(2*z2*w2^2);                     % determine the modal constants
Rm=A1./(w1^2-w.^2+j*2*z1*w*w1)+A2./(w2^2-w.^2+j*2*z2*w*w2); % estimated receptance from modal parameters 

%% plot the results
figure                                          % FRF modulus
plot(ff,20*log10(abs(Re)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(ff,20*log10(abs(Rm)),'--k','linewidth',3)
axis square,axis([0,40,-120,-20]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('|receptance| (dB ref 1 m/N)');

figure                                                 % phase
plot(ff,180/pi*unwrap(angle(Re)),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(ff,180/pi*unwrap(angle(Rm)),'--k','linewidth',4)
axis square,axis([0,40,-200,0]),grid
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('phase (degrees)');

figure                                       % Nyquist plot
plot(real(Re),imag(Re),'linewidth',4,'color',[0.6 0.6 0.6])
hold on
plot(real(Rm),imag(Rm),'--k','linewidth',4)
hold on
plot(real(R),imag(R),'k','linewidth',4)
axis square,axis([-0.01,0.01,-0.02,0]),grid
set(gca,'fontsize',16)
xlabel('real \{receptance\} (m/N)');
ylabel('imag \{receptance\} (m/N)');
