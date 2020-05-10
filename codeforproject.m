clear; clc; 
close all;
% Information given
h = 150;% w/m^2k convection coeff 
k = 11; % w/mk conduction coeff
N = 101; % nodes
r = 0.02; % meters
rhost = 7810; % density of steel
rhoair = 1.204; % density of air
V_room = 25; % meters cubed
cpst = 0.51; % kj/kgk specific heat steel
cpair = 1.006; % kj/kgk specific heat air
Tamb = 24; % room temperature centigrade
Tin = 350; % sphere temperature initial centigrade
Tfn = 24; % sphere temperature final centigrade
alpha = k/(rhost*cpst);
Ac = pi*r^2;
P = 4*r;
m=sqrt(h*P/(k*Ac));
% Creating grid necessary for nodal analysis
xi = linspace(0,r,N)';
dx = xi(2)-xi(1);

% Creating Matrix
A = zeros(N); 
A(1,1) = 1;


% Column Vector
b = zeros(N,1);
b(1) = -Tamb;


% NODES
for i = 2:N-1
    A(i,i-1) = alpha/dx^2;
    A(i,i) = -alpha/dx^2*(2+m^2*dx^2);
    A(i,i+1) = alpha/dx^2;
    b(i) = alpha*m^2*Tin;
end % This is for middle NODES

% Equations for next step
A(N,N-1) = 2*alpha/dx^2;
A(N,N) = -2*alpha/dx^2*(2+m^2*dx^2/2);
b(N) = alpha*m^2*Tin;

% Time variables and Steps
dt = 25e-3;
tf=60;
t = 0:dt:tf;
numTimehops = length(t);

% Matrices created
I = eye(N);
m1 = I-dt/2*A;
m2 = I+dt/2*A;

T = zeros(N,numTimehops);

T(:,1) = Tin;

invm1 = m1\I;

for i = 1:numTimehops-1
    T(1,i) = Tin;
    T(:,i+1) = invm1*m2*T(:,i)+dt*invm1*b;
    
end


%% Ploting
hold on;
grid on;
plot(xi*1000,T(:,1:100:end))
title(['T vs x step every', num2str(1000*dt),' from 0s to ', num2str(tf),'s']);
xlabel('x - [mm]'); 
ylabel('T - [^oC]');



