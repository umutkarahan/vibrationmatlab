clear all 

N=100;                                       % No. of samples per cycle
                                             % (choose 100, 2, 1.1, 0.9) 
df=2*pi/N;                                   % angular resolution
for phi=0:df:20*pi                           % angle
   X=exp(j*phi);                             % rotating vector
   x=[0 real(X)]; y=[0 imag(X)];             % x and y co-ordinates
   pause (1/(N))                             % creates a pause between plots
   
   th=0:0.01:2*pi;                           % create the circle
   Z=exp(j*th);
   zx=real(Z);zy=imag(Z);
   
   %% plot the results
   plot(zx,zy,'linewidth',3,'Color',[.6 .6 .6])
   hold on
   plot(x,y,'k','linewidth',4) 
   hold on
   plot(x,y,'ok','linewidth',2,'markerfacecolor','w','markersize',10)
   set(gca,'fontsize',16)
   axis([-1.2,1.2,-1.2,1.2])
   grid;axis square
   hold off
end