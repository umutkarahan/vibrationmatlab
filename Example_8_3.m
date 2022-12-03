clear all

%% parameters
m1=1;m2=1;m3=1;                               % masses
k1=1e4;k2=1e4;k3=1e4;k4=0*5e3;                % stiffness
M=[m1 0 0; 0 m2 0; 0 0 m3];                   % mass matrix
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];  % stiffness matrix

%% undamped natural frequencies
[V W]= eig (inv(M)*K);                        % calculation of eigenvalues
f123=sqrt(W)/(2*pi)                           % calculation of natural frequencies
%% anti-resonance frequencies X1/Fe1
M=[m2 0; 0 m3];                               % mass matrix
K=[k2+k3 -k3; -k3 k3];                        % stiffness matrix
[V1 W1]= eig (inv(M)*K);                      % calculation of eigenvalues                                        
a12=sqrt(W1)/(2*pi)                           % calculation of antiresonance frequencies

%% anti-resonance frequencies X2/Fe1
W2=(k3)/m1;                                   % calculation of the square of the natural frequency
a3=sqrt(W2)/(2*pi)                            % calculation of the anti-resonance frequency


