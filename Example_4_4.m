clear all 

%% Parameters
m = 1;                              % See Matlab example 3.1           
k=10000;                              
z = 0.1; c = 2*z*sqrt(m*k);          
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    
%% Time and frequency parameters
T=0.1;TT=1;                         % See Matlab example 3.1
fs=400;                             
df=1;                               
dt=1/fs; t=0:dt:T; tt=0:dt:TT;      
%% Theoretical IRF
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);   % IRF short duration windowd
hT=1/(m*wd)*exp(-z*wn*tt).*sin(wd*tt);% IRF long duration window

%% IRF plot
plot(tt,hT,'k','linewidth',4,'Color',[.6 .6 .6]) % plot IRF
hold on
plot(t,h,'k','linewidth',4)
grid; axis square 
axis([0,1,-1e-2 1e-2])
set(gca,'fontsize',24)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');

%% Theoretical FRF and sinc function
f=-fs/2:df:fs/2;                      % frequency vector                          % [rad/s]
w=2*pi*f;
H=1./(k-w.^2*m+j*w*c);                % FRF
fr=-1.5*fs/2:df:1.5*fs/2; 
W=T*sinc(T*fr).*exp(-j*pi*fr*T);      % sinc function
delta=sinc(fr).*exp(-j*pi*fr);        % approximate delta function

%% Convolve sinc function with FRF
C=conv(W,H)*df;Ca=conv(delta,H)*df;   % with and without truncation 
%% Plot results
fmin = 3*min(f); fmax = 3*max(f);     % set frequency range for animation     
Wf = fliplr(W);                       % flip the sinc function
ff = fliplr(-fr);
ff = ff + ( min(f)-max(ff) );         % slide range of W
fc = [ff f(2:end)];fc = fc+max(fr);   % range of convolved function
set(figure,'Position', [40, 40, 1450, 700]);  % set position of animation
  
subplot(2,1,1);
HN=abs(H)/max(abs(H));WN=abs(W)/max(abs(W));  % mormalised FRF and W
p = plot(f,HN,'k','linewidth',4); hold on     % plot of normalised FRF
gr=[.6 .6 .6];                                % define grey colour
q = plot(fr,WN,'k','linewidth',4,'Color',gr); % plot of normalised W
axis([fmin,fmax,0,1.1])                      
set(gca,'fontname', 'arial','fontsize',20)
xlabel('frequency (Hz)');
ylabel('normalised modulus');

% plot two vertical lines to show the range of ovelapped area
sl = line( [min(f) min(f)], [1.1 1.1], 'color','k');  % vertical line for overlap
hold on; grid on;
% ovelapped shaded region
sg = rectangle('Position', [min(f) 1 0 0], ...        % shaded region
                'FaceColor', [.9 .9 .9]);

subplot(2,1,2);
CdB=20*log10(abs(C));CadB=20*log10(abs(Ca));           % convolved values in dB
r = plot(fc,CdB,'linewidth',3,'Color',gr);hold on;     % plot of convolution
s = plot(fc,CadB,'k','linewidth',4);                   % as above no truncation
grid on; hold on;
axis([fmin,fmax,-140,max(CdB)+10])  
set(gca,'fontname', 'arial','fontsize',20)
xlabel('frequency (Hz)');
ylabel('modulus (dB ref 1N/m)');

%% animation block
for n=1:2/3*length(fc)                      
    pause(0);                                   % controls animation speed
    ff=ff+df;
    set(q,'XData',ff,'YData',WN);
        
    sx = min( max( ff(1), min(f) ), max(f) );   % left-hand boundary of overlap
    sx_a = [sx sx];
    set(sl,'XData',sx_a);
 
    ex = min( ff(end), max(f) );                % right-hand boundary of overlap
    ex_a = [ex ex];
    set(sl,'XData', ex_a);
       
    rpos = [sx 0 max(0.0001, ex-sx) 1.1];       % shading of overlap region
    set(sg,'Position',rpos);
    uistack(sg,'bottom');
      
    set(r,'XData',fc(1:n),'YData',CdB(1:n));    % plot of convolved function
    set(s,'XData',fc(1:n),'YData',CadB(1:n));   % as above no truncation
end;
