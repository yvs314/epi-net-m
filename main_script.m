% Networked SIR simulator

% Authorship: Shuang Gao, Rinel Foguen, and Yaroslav Salii
% This code depends on Trajectory.m for ploting Folder data for storing data 
% Folder /fig for saving the figures
% Folder /data is for loading and saving data
%
%% Clear the workspace
clc ; clear all; close all;
%% Set/Read the System Parameters

% Set the number of "non-terminal" compartments
% i.e., the number of all but the final R (recovered/removed) compartment
nodeDim = 2 ;  %number of states at each node (node's dimension)
         % first dimension is s, the fraction of susceptible
         % second dimension is i, the fraction of infected

% Set time horizon: [t_0, tFin]
% t_0 = 0, %time starts at 0
tFin = 2;          %end time (finite horizon); say, 30 days

%Time Discretization Parameters
%TODO: move on to steps per day
nSteps = 1000 ;       %number of time steps
dt = tFin/nSteps ;      % time step size (unifrom)

%Model Parameters
gamma = 1/8.3 ; %removed rate at each node
beta = 1/2.5 ; %infectious rate at each node

%% Read the Initial Values (inc. node number)
%a .CSV with the cols {AP_ID,AP_code,N_i,S_i,I_i,R_i,City_name},
%each row defines a node
tInitialVals = readtable("data/init-20-v3.csv");
nodeNum = size(tInitialVals,1); %as many nodes as there are rows
%make a population vector bN out of tIVs' 3rd column N_i
bN = table2array(tInitialVals(:,3));
%read the absolute S_i and I_i
X0_absolute = table2array(tInitialVals(:,[4,5]));
%now make a fractional version: s_i = S_i / N_i, etc.
X0_frac = X0_absolute ./ bN;

%re-shape X0_frac from [n\times2] into column vector [2n\times1]
%of the form [s1;i1;s2;i2...s_n;i_n]
X_0 = flattenRowMjr(X0_frac);

%% Read the Coupling Data (daily passengers)
DATA = load("data/10000000-AP-daily-flights-BTS-2019.dat");
%let us first test the no-travel case
A = zeros(nodeNum,nodeNum);
%A = ones(m,m);
%A = DATA/norm(DATA, 'inf')*100; 
v_forDiag = ones(nodeNum,1); %[m\times 1], column vector
%make sure A's diagonal elements are mostly 1
%YS: note that I HAET this idea, I'd rather have a separate term for a pop's own attack rate
A = purgeDiag(A)+diag(v_forDiag);

%% Set The Control Parameters
u = zeros(1,nodeNum) ; %load the control here

%% Define the Dynamics, set the Initial Values
% initialize the node states' time series X with zeros at all time steps
X = zeros(nodeDim*nodeNum,nSteps+1) ;
% set the *initial conditions*
X(:,1) = X_0; 

%Build linear drift term YS:What?
D1 = zeros(nodeDim,nodeDim) ; 
D1(2,2) = - gamma ;

%Build the nonlinear drift term YS:What?
D2 = zeros(nodeDim,nodeDim) ;
D2(1,2) = -1 ;
D2(2,2) = 1 ;

%% Run and Time the Dynamic Simulation
tic
for j=1:nSteps %for every time step
    
    % Build the diagonal matrix with suspected states only
    v1 = zeros(1,nodeNum) ;
    for k=1:2:(nodeDim*nodeNum) %What?
        v1 = X(k,j) ;
    end
    D3 = diag(v1) ;
    U = diag(u) ;
    
    %Compute the time series
    %Z =  kron(eye(m),D1).*dt +  dt.*beta.*kron( (D3.*A) , D2 ) ...
    %     -  dt.*beta.*kron( (U.*D3.*A) , D2 ) ;
     
    Z =  kron(eye(nodeNum),D1)*dt +  dt*beta*kron( (D3*A) , D2 ) ...
         -  dt*beta*kron( (U*D3*A) , D2 ) ;  
     % .* is element-wise multiplication 
     % we are using matrix multiplication here
    
    for h=1:(nodeDim*nodeNum)
        X(h,j+1) = X(h,j) +  Z(h,:)*X(:,j) ;
    end
    
end
toc

%% Analysis and output
% prep the time for graphing
tSpan = [0,tFin];
tSteps = nSteps+1;
t      = linspace(tSpan(1),tSpan(2),tSteps);

nAgent =  nodeNum;
nState =  nodeDim;
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
