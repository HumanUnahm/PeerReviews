clear; clc; 
close all;
% Information given
h = 150;% w/m^2k convection coeff 
k = 11; % w/mk conduction coeff
N = 101; % nodes
r = 0.2/2; % meters % Bad Conversion cm to meters % Used diameter (0.2m) instead of radius (0.1m)
rhost = 7810; % density of steel % Bad Conversion (kg is SI unit. grams is not standard SI).
rhoair = 1.204; % density of air
V_room = 25; % meters cubed
cpst = 0.51*1e3; % kj/kgk specific heat steel % Bad Conversion (J is SI unit, kJ is not standard SI)
cpair = 1.006; % kj/kgk specific heat air
Tamb = 24; % room temperature centigrade
Tin = 350; % sphere temperature initial centigrade
Tfn = 24; % sphere temperature final centigrade
alpha = k/(rhost*cpst);
Ac = pi*r^2;
P = 4*r;
m=sqrt(h*P/(k*Ac));
% Creating grid necessary for nodal analysis
ri = linspace(0,r,N)';
dr = ri(2)-ri(1);

% Creating Matrix
A = zeros(N); 
A(1,1) = -6*alpha/dr^2;
A(1,2) = 6*alpha/dr^2;


% Column Vector
 b = zeros(N,1);
% b(1) = -Tamb;


% NODES
for i = 2:N-1
    A(i,i-1) = alpha*(1-dr/ri(i))/dr^2;
    A(i,i) = -2*alpha/dr^2;
    A(i,i+1) = alpha*(1+dr/ri(i))/dr^2;
    %b(i) = alpha*m^2*Tin;
end % This is for middle NODES


% Equations for next step
A(N,N-1) = 2*alpha/dr^2;
A(N,N) = -2*alpha*(1+(1+dr/r)*h*dr/k)/dr^2;
b(N) = -2*alpha*(1+dr/r)/dr^2;


% Time variables and Steps
dt = 25e-3;
tf=300;
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
    T(:,i+1) = invm1*m2*T(:,i)+dt*invm1*b;
    
end


%% Plotting 0 to 60s
hold on;
grid on;

timechecks = [1,2,5,10,20,30,45,60]*(1/dt); % Replicates Kevin's Results Gifts for temperature profiles up to one minute. 
timechecks(1) = 1;
plot(ri*1000,T(:,timechecks)); % m to mm conversion is 1000
%title(['T vs x step every', num2str(100*dt),' from 0s to ', num2str(tf),'s']);
title('replicating Kevin Gifted results');
xlabel('x - [mm]'); 
ylabel('T - [^oC]');
legend

%% Plotting 0 to 1 hour

% Time variables and Steps
%dt = 25e-3;
%tf2=3600;
%t = 0:dt:tf;
%numTimehops = length(t);

% Matrices created
%I = eye(N);
%m1 = I-dt/2*A;
%m2 = I+dt/2*A;

%T = zeros(N,numTimehops);

%T(:,1) = Tin;

%invm1 = m1\I;

%for i = 1:numTimehops-1
 %   T(:,i+1) = invm1*m2*T(:,i)+dt*invm1*b;
    
%end


%hold on;
%grid on;

%% Jacobi Solution for equation
%tic
%Tsolg = Gaus(A,b)
%Guasselim_time = toc

% Jacobi time
%tol = 0.125e-12;
%itmax = 1e6;
%Tguess = 0*ones(N,1);
%tic

%Tsol_cobi = Jacobmethod(A,b,Tguess,tol,itmax);
%Jacobi_time = toc






