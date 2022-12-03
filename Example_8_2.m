clear all

%% parameters
m1=1;m2=1;m3=1;                               % masses
k1=1e4;k2=1e4;k3=1e4;k4=0*5e3;                % stiffness
M=[m1 0 0; 0 m2 0; 0 0 m3];                   % mass matrix
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];  % stiffness matrix

%% undamped natural frequencies
[V W]= eig (inv(M)*K);                         % calculation of eigenvalues and eigenvectors
R=sqrt(W)/(2*pi)                               % Calculation of natural frequencies
V1=V(:,1)/max(abs(V(:,1)))                     % 1st mode shape
V2=V(:,2)/max(abs(V(:,2)))                     % 2nd mode shape
V3=V(:,3)/max(abs(V(:,3)))                     % 3rd mode shape

