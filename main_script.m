%% Networked SIR simulator
%
% Authorship: Shuang Gao, Rinel Foguen, and Yaroslav Salii
% 
% Plotting facilities: 
    %Trajectory.m (SG) for 3D overview
    %figStacked.m (YS) for per-node stacked i+s+r 
% Directories:
    %./fig is for saving the figures (not tracked by git)
    %./data is for initial values etc. (instance collection)
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
tFin = 150;          %end time (finite horizon); say, 30 days

%Time Discretization Parameters
%TODO: move on to steps per day
stepsPerDay=1;
nSteps = tFin * stepsPerDay;       %number of time steps
dt = tFin/nSteps ;      % time step size (unifrom)

%Model Parameters
gamma = 1/8.3 ; %removed rate at each node
beta = 1/2.5 ; %infectious rate at each node

%% Read the Initial Values (inc. node number)
%a .CSV with the cols {AP_ID,AP_code,N_i,S_i,I_i,R_i,City_name},
%each row defines a node
iValPath="data/init-4-first.csv"; %where do we keep the IVs,
%iFlugPath="data/flug-2-first.dat"; %where do we keep daily passengers

tInitialVals = readtable(iValPath);
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

%% Set the Coupling Data (daily passengers) / dbg: set to no-travel or eps-one
%DATA = load(iFlugPath); %daily passengers, from a file
%A = eye(nodeNum); %no connection
A = mkEpsOneMx(nodeNum,1e-6); % epsilon-one full connectivity

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

% Correct implementation
    % pick S
    v1  = kron(eye(nodeNum),[1,0])* X(:,j);
    D3  = diag(v1);
    U   = diag(u) ;  
    Z   =  kron(eye(nodeNum),D1)*dt +  dt*beta*kron( (D3*A) , D2 ) ... 
        -  dt*beta*kron( (U*D3*A) , D2 ) ;
    X(:,j+1) = X(:,j) +  Z*X(:,j);
end
toc

%% Analysis and output
% prep the time steps as a series (row) for graphing
tSpan = [0,tFin];
tSteps = nSteps+1;
t      = linspace(tSpan(1),tSpan(2),tSteps);

%original X is [nodeDim \cdot nodeNum \times tSteps],
%[ s_1(0) s_1(1) ... s_1(t_f) 
% ; i_1(0) i_1(1) ... i_1(t_f)
% ; s_2(0) s_2(1) ... s_2(t_f) 
% ; i_2(0) i_2(1) ... i_2(t_f) ]

% [s,i,r]Evo are [nodeNum \times tSteps], e.g., sEvo[k,t] is s_k(t)
%carve the s and i compartments' time series out of the whole X
sEvo = X(1:nodeDim:nodeNum*nodeDim,:);
iEvo = X(2:nodeDim:nodeNum*nodeDim,:);
% r is not computed initially, so we add it by "closure"
rEvo = ones(size(sEvo)) - sEvo - iEvo;

%poor, hacky, self-repeating construction of absolute values
S_Evo = sEvo .* bN;
I_Evo = iEvo .* bN;
R_Evo = rEvo .* bN;
%Figures are meant to be  written into ./figDirName; make sure it exists
figDirName = 'fig';
if(~exist(figDirName,'dir'))
    mkdir(figDirName);
end


f =  tiledlayout('flow'); %built-in layout, aimes at 4:3 for the tiles
for thisNode = 1:nodeNum
    nexttile
    figStacked(t,sEvo(thisNode,:) ...
                       ,iEvo(thisNode,:) ...
                       ,rEvo(thisNode,:) ...
                       ,tInitialVals.City_name(thisNode));
end

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

%% aux: Make Epsilon-One Coupling Matrix
% primitive coupling: \epsilon\in[0,1] off-diagonal, 1 on diagonal
function Aout = mkEpsOneMx(size, eps)
    Aout= purgeDiag(repmat(eps,size))+eye(size);   
end