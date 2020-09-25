%% Networked SIR simulator
%
% Authorship: Shuang Gao, Rinel Foguen, and Yaroslav Salii
% 
% Plotting facilities: 
    %Trajectory.m (SG) for 3D overview
    %figStacked.m (YS) for per-node stacked i+s+r 
% Directories: (set by $figDir and $instDir)
    %./fig is for saving the figures (not tracked by git!)
    %./data is for initial values etc. (instance collection)
% Naming conventions: fields separated by "-" subfields by "_" 
% Instances:    (with string $NAME and integer $SIZE)
    % $NAME_$SIZE-init.csv holds initial conditions; 
    % $NAME_$SIZE-flug.dat holds daily passengers matrix
    %sample: first4-init.csv IVs, Atlanta, Boston, Charlotte, and Denver
            %first4-flug.dat daily psg., Atlanta, Boston, Charlotte, and Denver
%  Output figures: ($CPG: [flug | eps-one] coupling type, $TYPE:[OVR|STK] figure type)
    %$NAME_$SIZE-$CPG-$TYPE.$EXT
    %OVR: $NAME_$SIZE-$CPG-OVR.png All-nodes Overview (via Trajectory.m)
    %STK: $NAME_$SIZE-$CPG-STK.png Stacked per-node i+s+r, all on one figure (via figStacked.m)
%  Output tables: $NAME_$SIZE-$CPG-[abs | frac]-out.csv

%consider adding $beta, $gamma, $tFinal, and $stepsPerDay
%% Clear the workspace
clc ; clear all; close all;
%% Set/Read the System Parameters

%set the instance location
%instDir="data/by-tract"; %by-tract are FluTE-derived, larger instances
instDir="data"; %just the small airport-derived cases here

%initial values instances have the form $instName-init.csv

%set the instance name
%instName="seattle_124"; %smalles of the FluTe-derived
instName="first_4"; %toy 4-city AP-derived instance, Atlanta, Boston, Charlotteville, Denver


% if false, use eps_one coupling
% if true, read the daily passengers table from $instName-flug.dat
useFlightData=false;

% Set the number of "non-terminal" compartments
% i.e., the number of all but the final R (recovered/removed) compartment
nodeDim = 2 ;  %number of states at each node (node's dimension)
         % first dimension is s, the fraction of susceptible
         % second dimension is i, the fraction of infected

% Set time horizon: [t_0, tFin]
% t_0 = 0, %time starts at 0
tFin = 150;          %end time, in days (finite horizon);

%Time Discretization Parameters
stepsPerDay=2; %must be at least 1; better be integer
nSteps = tFin * stepsPerDay;       %number of time steps
dt = tFin/nSteps ;      % time step size (uniform)

%Model Parameters
gamma = 1/8.3 ; %removed rate at each node
beta = 1/2.5 ; %infectious rate at each node

%% Set the IO parameters: input instances
fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields
IV_suff="init.csv"; %all instances IVs end like this 
flug_suff="flug.dat"; %all instances daily psgrs end like this

%TODO: add error if IVs not found
iValPath=fullfile(instDir,instName+fnSep+IV_suff); %path to the IVs
iFlugPath=fullfile(instDir,instName+fnSep+flug_suff); %path to the flight data, if any

%% Set the IO parameters: output location and naming
outDir = "out"; %write the output figures and tables here
if(~exist(outDir,"dir"))
    mkdir(outDir); %make sure it exists
end

% set the right coupling code
nCPG="undef"; %init the coupling code; I'm scared of var init inside if block
if (useFlightData)
    nCPG = "flug";
else
    nCPG="eps_one";
end
%nCPG=(useFlightData)?("flug"):("eps_one"); %default to epsilon-one coupling
%brief: $NAME-$CPG-[STK|OVR].$ext STK for Stacked, OVR for Overview
pathOutFigStacked = fullfile(outDir,join([instName,nCPG,"STK"],fnSep));
pathOutFigOverview = fullfile(outDir,join([instName,nCPG,"OVR"],fnSep));

pathOutTabAbs = fullfile(outDir,join([instName,nCPG,"abs"],fnSep));
pathOutTabFrac = fullfile(outDir,join([instName,nCPG,"frac"],fnSep));
%$NAME_$SIZE-$CPG-beta_$beta-gamma_$gamma-T_$tFinal.$EXT

%% Read the Initial Values (inc. node number)

%a .CSV with the cols {AP_ID,AP_code,N_i,S_i,I_i,R_i,City_name},
%each row defines a node
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

%% Read/Set the Coupling Data
%CAVEAT: no symmetrization yet
if(useFlightData)  %read and process daily passengers
    % must have city pops N=(N1,...,Nn) on the diagonal
    A = purgeDiag(load(iFlugPath))+diag(bN);
else
    A = mkEpsOneMx(nodeNum,1e-4); % epsilon-one full connectivity
   %A = eye(nodeNum); %no coupling, another debug variant
end

%% Set The Control Parameters
u = zeros(1,nodeNum) ; %load the control here

%% Define the Dynamics, set the Initial Values
% initialize the node states' time series X with zeros at all time steps
X = zeros(nodeDim*nodeNum,nSteps+1) ;
% set the *initial conditions*
X(:,1) = X_0; 

%Build the Infected -> Removed term (the *linear* drift term)  
% D1=[0 0;0 \gamma]
D1 = zeros(nodeDim,nodeDim); 
D1(2,2) = - gamma ;

%Build the Susceptible -> Infected term (the *nonlinear* drift term)
% D2 = [0 -1; 0 1]
D2 = zeros(nodeDim,nodeDim) ;
D2(1,2) = -1 ;
D2(2,2) = 1 ;

%if using (absolute) daily passengers
if(useFlightData)
    %add division by each city's population
    absCorr=diag(arrayfun(@(x) 1/x,bN));
else %we use fractional "connection intensities", not absolute psg/day
    %no need to divide by anything, just use the scalar 1 matrix
    absCorr=eye(nodeNum);
end

%% Run and Time the Dynamic Simulation
disp("Running the simulation");
tic
for j=1:nSteps %for every time step
    %form the infection propagation and transport term D3:
    %from X=[s_1(j) i_1(j) ... s_last(j) i_last(j)],
    %kron() picks every second one, the susceptibles [s_1(j) ... s_last(j)]
    %if using *absolute passengers per day* (useFlightData=true)
    %then we also divide each by the city population (see absCorr above)
    D3  = diag(kron(absCorr,[1,0])* X(:,j)); 
    U   = diag(u) ;  
    Z   =  kron(eye(nodeNum),D1)*dt +  dt*beta*kron( (D3*A) , D2 ) ... 
        -  dt*beta*kron( (U*D3*A) , D2 ) ;
    X(:,j+1) = X(:,j) +  Z*X(:,j);
end
toc
disp("Numerics complete");
%% Analysis and output
disp("Prep for graphing");
tic
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
toc
%% Plot the 2D Stacked Per-Node i+s+r (all fractional)
disp("Plotting the 2D stacked per-node i+s+r");
tic
%only so many fit on a tiled figure; 20 ok on my 13" display
%but they can't reasonably be printed to a single .png, print at most 9
maxTiles = 20;
% make a figure object for these tiled per-node stacked plots
fStacked = figure("name","Stacked Per-Node i+s+r, Fractional");
% start tiling
tiledlayout("flow"); %built-in layout, aimes at 4:3 for the tiles
for thisNode = 1:min(maxTiles,nodeNum)
    nexttile
    figStacked(t,sEvo(thisNode,:) ...
                       ,iEvo(thisNode,:) ...
                       ,rEvo(thisNode,:) ...
                       ,tInitialVals.(7)(thisNode));
end
toc
%% Plot the 3D absolute evolution
disp("Plotting All-nodes S,I,R Trajectories");
tic
% make a figure object for the tiled abs-all-nodes plots
fAllNodesAbs = figure("name","S, I, and R for All Nodes Together, Absolute Values.");
%prep the nodes/APs labels for the plots
nodeLabels = table2array(tInitialVals(:,7))';
%start tiling
tiledlayout("flow");
nexttile;
fS=Trajectory(S_Evo,t,"b","abs-Susceptible",nodeLabels);
nexttile;
fI=Trajectory(I_Evo,t,"r","abs-Infected",nodeLabels);
nexttile;
fR=Trajectory(R_Evo,t,"m","abs-Removed",nodeLabels);
nexttile;
%total population: debug value only; assert constant populations
f3All = Trajectory(S_Evo+I_Evo+R_Evo,t,"k","abs-Total",nodeLabels);
toc
%% Figure export, default to png
%print the big picture
print(fAllNodesAbs,pathOutFigOverview,'-dpng');
%print the tiled stacked per-nodes (just the first maxTiles)
print(fStacked,pathOutFigStacked,'-dpng');
%% Time series export, to a .csv
disp("Re-shaping the output time series");
tic
%generate the all-ticks output tables
outTabFrac=mkOutTable(sEvo,iEvo,rEvo);
outTabAbs=mkOutTable(S_Evo,I_Evo,R_Evo);
toc

tic
%generate the per-day output tables 
outTabFracDaily=mkOutTable(sEvo(:,1:stepsPerDay:end) ...
                        ,iEvo(:,1:stepsPerDay:end) ...
                        ,rEvo(:,1:stepsPerDay:end));

outTabAbsDaily=mkOutTable(S_Evo(:,1:stepsPerDay:end) ...
                        ,I_Evo(:,1:stepsPerDay:end) ...
                        ,R_Evo(:,1:stepsPerDay:end));
toc
%write the output tables to .csv
writematrix(outTabFracDaily,strcat(pathOutTabFrac,'.csv'))
writematrix(outTabAbsDaily,strcat(pathOutTabAbs,'.csv'))

%% dbf lf: for calling local functions from command line
lf = localfunctions;
%% aux: Make Output Table
%TODO: switch the output display from ticks to days, using stepsPerDay
% takes the 3 compartments' time series, e.g. s_1,\dots,s_n[0,tFin], i, r
% and puts them into a table of the form: \forall t
% t, s_1(t), i_1(t), r_1(t), ... , s_n(t), i_1(t), r_n(t)
% NB!: assume time series are of the size [nodeNum\times nSteps]
function outT = mkOutTable(S, I, R)
nNodes = size(S,1);    %must match for S, I, and R
nTics = size(S,2); %must match for S, I, and R

tmpA=[flattenRowMjr(S)' flattenRowMjr(I)' flattenRowMjr(R)' ];

outT = zeros(nTics,3*nNodes +1);  
outT(:,1)=0:(nTics-1); %write time ticks as the first column, from tick 0
for k = 1:nNodes %now copy the info from tmpA with a sliding window
    outT(1:nTics,3*(k-1)+2:3*k+1) =  tmpA( (k-1)*nTics +1:k*nTics,:);
end
end %end function, return outT
%% aux: Flatten Row-Major,
% turns [N\times nCol] matrix into a column vector [nCol*N\times 1]
% in row-major order, e.g. [s1 i1; s2 i2] -> [s1; i1; s2; i2]
function Xout = flattenRowMjr(Xin)
    Xout = reshape(Xin',1,numel(Xin));
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