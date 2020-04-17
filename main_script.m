% Authorship: Shuang Gao and Rinel Foguen
% vandalized by Yaroslav Salii
% This code depends on Trajectory.m for ploting Folder data for storing data 
% Folder /fig for saving the figures
% Folder /data is for loading and saving data

%
%% Clear the workspace
clc ; clear all; close all;
%% System Parameters

tic
%Time Discretization Parameters
%TODO: move on to steps per day
T = 60;          %continuous time duration (30 days)
N = 100 ;       %number of time steps
dt = T/N ;      % time step size (unifrom)

%Model Parameters
m = 20 ; % number of nodes in the network
d = 2 ;  %dimension of the states at each node
         % first dimension is number of susceptible
         % second dimension is number of infected

gamma = 1/8.3 ; %removed rate at each node
beta = 1/2.5 ; %infectious rate at each node
X = zeros(d*m,N+1) ; %states time series 

%% Data input block: read the flight data and the initial conditions
%read the inital values & info on them,
%a .CSV with the cols {AP_ID,AP_code,N_i,S_i,I_i,R_i,City_name}
tInitialValInfo = readtable("data/init-20-v3.csv");
%make a population vector bN out of its 3rd column N_i
bN = table2array(tInitialValInfo(:,3));
%read the absolute S_i and I_i
X0_absolute = table2array(tInitialValInfo(:,[4,5]));
%now make a fractional version: s_i = S_i / I_i, etc.
X0_frac = X0_absolute ./ bN;

%re-shape X0_frac from [20x2] into column vector [40x1]
%of the form [s1;i1;s2;i2...sN;iN]
X_0 = flattenRowMjr(X0_frac);

%% just read the flight data
DATA = load("data/10000000-AP-daily-flights-BTS-2019.dat");

%% System Data
%let us first test the no-travel case
A = zeros(m,m);
%A = ones(m,m);
%A = DATA/norm(DATA, 'inf')*100; 
v_forDiag = ones(m,1); %[m\times 1], column vector
%make sure A's diagonal elements are mostly 1
%YS: note that I HAET this idea, I'd rather have a separate term for a pop's own attack rate
A = purgeDiag(A)+diag(v_forDiag);

%% Control Parameters

u = zeros(1,m) ; %load the control here
%% Dynamic Simulation

%initializing the time series
X(:,1) = X_0; 


%Build linear drift term YS:What?
D1 = zeros(d,d) ; 
D1(2,2) = - gamma ;

%Build the nonlinear drift term YS:What?
D2 = zeros(d,d) ;
D2(1,2) = -1 ;
D2(2,2) = 1 ;

for j=1:N %for every time step
    
%     Build the diagonal matrix with suspected states only
%     v1 = zeros(1,m) ;
%    
%     for k=1:2:(d*m) %What?
%         v1 = X(k,j) ;
%     end
% 
%     D3 = diag(v1) ;
%     U = diag(u) ;
%     % 
%     
% %%YS: fractional input assumed (\beta_i is not divided by pop_i)
%     %Compute the time series
%     Z =  kron(eye(m),D1).*dt +  dt.*beta.*kron( (D3.*A) , D2 ) ...
%          -  dt.*beta.*kron( (U.*D3.*A) , D2 ) ;
%        
%     for h=1:(d*m)
%         X(h,j+1) = X(h,j) +  Z(h,:)*X(:,j) ;
%     end

% Correct implementation
    % pick S
    v1  = kron(eye(m),[1,0])* X(:,j);
    D3  = diag(v1);
    U   = diag(u) ;
    
    Z   =  kron(eye(m),D1)*dt +  dt*beta*kron( (D3*A) , D2 ) -  dt*beta*kron( (U*D3*A) , D2 ) ;
    X(:,j+1) = X(:,j) +  Z*X(:,j);
end
toc
tSpan = [0,T];
tSteps = N+1;
t      = linspace(tSpan(1),tSpan(2),tSteps);

%% Analysis and output
nAgent =  m;
nState =  d;
sEvo = X(1:2:nAgent*nState,:);
iEvo = X(2:2:nAgent*nState,:);
rEvo = ones(size(sEvo)) - sEvo - iEvo;

%Figures are written into ./figDirName
figDirName = 'fig';
if(~exist(figDirName,'dir'))
    mkdir(figDirName);
end
%TODO: refit Trajectory calls to use figDirName
Trajectory(sEvo,t,'b','SIR-Susceptible','fig/SIR-Susceptible')
Trajectory(iEvo,t,'r','SIR-Infected','fig/SIR-Infected')
Trajectory(rEvo,t,'m','SIR-Removed','fig/SIR-Removed')
Trajectory(rEvo+iEvo+sEvo,t,'k','Total','fig/SIR-Total')


%% aux: Flatten Row-Major,
% turns [N\times nCol] matrix into a column vector [nCol*N\times 1]
% in row-major order, e.g. [s1 i1; s2 i2] -> [s1; i1; s2; i2]
function Xout = flattenRowMjr(Xin)
   tXin = Xin';
   Xout = tXin([1:numel(tXin)]);
end

%% aux: Purge Diagonal
% make sure the diagonal is all zeros
function Aout = purgeDiag(Ain)
    Aout = Ain - diag(diag(Ain));
end
