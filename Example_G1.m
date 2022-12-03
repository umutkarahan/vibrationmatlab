clear all 

%% parameters
m = 1;                              % [kg]       see MATLAB Example 3.1          
k=10000;                            % [N/m]  
z = 0.1; c = 2*z*sqrt(m*k);         % [Ns/m] 
wn=sqrt(k/m); wd=sqrt(1-z^2)*wn;    % [rad/s]

%% time and frequency parameters
fs=500;                             % [Hz]      see MATLAB Example 3.1 
T=0.5; dt=1/fs; t=0:dt:T;           % [s]

%% IRF
ho=1/(m*wd)*exp(-z*wn*t).*sin(wd*t);           % IRF
hn=ho/max(ho);                                 % normalised IRF

tmin = min(t)-abs(max(t)-min(t))-0.2;          % time range for the plot
tmax = max(t)+abs(max(t)-min(t))+0.2;

N=1;N=200;                                     % use N=1 for s scled delta function and N=200 for step changes in the force input
N1=round(length(t)/10); N2=length(t)-N1-N;
fo=[zeros(1,N1) ones(1,N) zeros(1,N2)];        % force input 

%% convolution
c = dt*conv(ho,fo);                            % convolution

%% animation
h = fliplr(hn);                                % flip the IRF
tf = fliplr(-t);
tf = tf + ( min(t)-max(tf) );                  % time window

tc = [ tf t(2:end)];                           % time range of input
tc = tc+max(t);
set(figure,'Position', [40, 40, 1450, 700]);   % figure position
gr=[.6 .6 .6];                                 % grey shade

%% plot the results
ax = subplot(2,1,1);                           % first plot
p = plot(tf, h,'k','linewidth',3); hold on
q = plot(t, fo,'linewidth',3,'Color',0.6*gr); 
axis([tmin,tmax,1.2*min(hn),1.2*max(hn)])  
ym = get(ax, 'ylim');
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('normalised IRF and force');
%%
sl = line( [min(t) min(t)], [ym(1) ym(2)], 'color','k');     % vertical lines for shaed region 
el = line( [min(t) min(t)], [ym(1) ym(2)], 'color', 'k');
hold on; grid on;
sg = rectangle('Position', [min(t) ym(1) 0 ym(2)-ym(1)], ... % shaded region
                 'FaceColor', [.9 .9 .9]);

ax2 = subplot(2,1,2);                                        % second plot
r = plot(tc,c,'k','linewidth',3);grid on; hold on;
s = plot(tc,c,'linewidth',3,'Color',gr);grid on; hold on;
uistack(s,'bottom');
axis([tmin,tmax,1.2*min(c),1.2*max(c)])  
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement (m)');

%% animation
for n=1:16/16*length(tc)                          
    pause(0.01);                                  % controls animation speed
    tf=tf+dt;
    set(p,'XData',tf,'YData',h);
    sx = min( max( tf(1), min(t) ), max(t) );     % left hand boundary of overlap
    sx_a = [sx sx];
    set(sl,'XData',sx_a);                         % right hand boundary of overlap
    ex = min( tf(end), max(t) );  
    ex_a = [ex ex];set(el,'XData', ex_a);         % shading of overlap region
    rpos = [sx ym(1) max(0.0001, ex-sx) ym(2)-ym(1)];  
    set(sg,'Position',rpos);
    uistack(sg,'bottom');
    set(r,'XData',tc(1:n),'YData',c(1:n));        % plot convolved function
end