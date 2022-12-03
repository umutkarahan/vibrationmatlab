clear all 

%% parameters
m = 1;                              % [kg]        see Matlab Example 3.1   
k=10000;                            % [N/m]  
z = 0.1; c = 2*z*sqrt(m*k);         % [Ns/m] 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % [rad/s]
fd=wd/(2*pi);
%% FRFs
df=0.4;f=-100:df:100;                           % [Hz]      frequency vector
w=2*pi*f;                                       % [rad/s]
H=1./(z*wn+j*w);                                % [m/N]     FRF of envelope
HH=1./(k-w.^2*m+j*w*c);                         % FRF of SDOF system
fr=-0.2*100:df:0.2*100;                         % frequency vector
W=zeros(length(fr),1);A=1/(2*j)*1/(m*wd);       % W is the FRF of the oscilatory term. It has units of m/N
W((length(fr)+1)/2-round(fd/df))=1/df*A;
W((length(fr)+1)/2+round(fd/df))=-1/df*A;

%% convolution
C=conv(W,H)*df;

%% plot the results 
fmin = 1.5*min(f); fmax = 1.5*max(f);            % set the frequency range for the graph
Wf = fliplr(W);                                  % flip the envelope
ff = fliplr(-fr);
ff = ff + ( min(f)-max(ff) );                    % slide range of W
fc = [ff f(2:end)];fc = fc+max(fr);              % range of convolved function
set(figure,'Position', [40, 40, 1450, 700]);     % set position of the animation
  
subplot(2,1,1);                                  % plot of separate functions
HN=abs(H)/max(abs(H));WN=abs(W)/max(abs(W));     % normalised functions
p = plot(f,HN,'k','linewidth',3); hold on        % normnalised FT of the envelope
gr=[.6 .6 .6];                                   % define grey colour
q = plot(fr,WN,'k','linewidth',3,'Color',0.6*gr); % normalised W 
axis([fmin,fmax,0,1.1])  
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('normalised modulus');
sl = line( [min(f) min(f)], [1.1 1.1], 'color','k'); % vertical line for overlap
hold on; grid on;
sg = rectangle('Position', [min(f) 1 0 0], ...       % shaded region
                'FaceColor', [.9 .9 .9]);
            
subplot(2,1,2);                                      % plot of convolved values in dB
CdB=20*log10(abs(C));HHdB=20*log10(abs(HH));
r = plot(fc,CdB,'k','linewidth',3);hold on
p = plot(f,HHdB,'k','linewidth',3,'Color',gr); 
hold on;grid on;
axis([fmin,fmax,-120,max(CdB)+10])  
set(gca,'fontsize',16)
xlabel('frequency (Hz)');
ylabel('modulus (dB ref 1N/m)');

%% animation block
for n=1:length(fc)
    pause(0);                                   % controls animation speed
    ff=ff+df;
    set(q,'XData',ff,'YData',WN);               
    sx = min(max(ff(1),min(f)),max(f));         % left-hand boundary of overlap  
    sx_a = [sx sx];
    set(sl,'XData',sx_a);
    ex = min( ff(end), max(f) );                % right hand boundary of overlap
    ex_a = [ex ex];
    set(sl,'XData', ex_a);
    rpos = [sx 0 max(0.0001, ex-sx) 1.1];       % shading of overlap region
    set(sg,'Position',rpos);
    uistack(sg,'bottom');
    set(p,'XData',f,'YData',HHdB);              % plot of FRF
    uistack(p,'bottom');
    set(r,'XData',fc(1:n),'YData',CdB(1:n));    % plot of convolved function
end