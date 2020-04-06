%% Authorship ...
% This code depends on Trajectory.m for ploting
% Folder data for storing data
% Folder fig for saving the figures

%

%% Clear the workspace 
clc ;
clear all
close all

%% System Parameters

tic
%Time Discretization Parameters
T = 2;          %continuous time duration
N = 1000 ;       %number of time steps
dt = T/N ;      % time step size (unifrom)

%Model Parameters
m = 20 ; % number of nodes in the network
d = 2 ;  %dimension of the states at each node
         % first dimension is number of suspected
         % second dimension is number of infected

gamma = 1 ; %removed rate at each node
beta = 1 ; %infectious rate at each node
X = zeros(d*m,N+1) ; %states time series 

% Load network data
DATA = load('data/10000000-AP-flights-BoTS-2019.dat');

%% System Data

%A = zeros(m,m); %load the contact network here
%A = ones(m,m);
A =DATA/norm(DATA, 'inf')*10; 
X_0 = ones(d*m,1) ; %load initial states at all nodes
%X_0 = randperm(,)
%X_0 = [10,0,10,2,10,1,10]' ;

%% Control Parameters
u = zeros(1,m) ; %load the control here

%% Dynamic Simulation

%initializing the time series
X(:,1) = X_0; 

%Build linear drift term
D1 = zeros(d,d) ; 
D1(2,2) = - gamma ;

%Build the nonlinear drift term
D2 = zeros(d,d) ;
D2(1,2) = -1 ;
D2(2,2) = 1 ;

for j=1:N
    
    % Build the diagonal matrix with suspected states only
    v1 = zeros(1,m) ;
    for k=1:2:(d*m)
        v1 = X(k,j) ;
    end
    D3 = diag(v1) ;
    U = diag(u) ;
    
    
    %Compute the time series
    Z =  kron(eye(m),D1).*dt +  dt.*beta.*kron( (D3.*A) , D2 ) -  dt.*beta.*kron( (U.*D3.*A) , D2 ) ;
    
    for h=1:(d*m)
        X(h,j+1) = X(h,j) +  Z(h,:)*X(:,j) ;
    end
    
end
toc
tSpan = [0,T];
tSteps = N+1;
t      = linspace(tSpan(1),tSpan(2),tSteps);

nAgent =  m;
nState =  d;
sEvo = X(1:2:nAgent*nState,:);
IEvo = X(2:2:nAgent*nState,:);
Trajectory(sEvo,t,'SIR-Sucieptible','fig/SIR-Sucieptible')
Trajectory(IEvo,t,'SIR-Infected','fig/SIR-Infected')
