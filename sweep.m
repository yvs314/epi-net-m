%% Forward-Backward Sweep for Net-SIR with Social Distancing Control

%Author: Yaroslav Salii, 2021+
% Directories: 
    %./data/by-tract is for problem instances (IVs + travel matrix)
    %./fig is for the figures (not tracked by git)
    %./out is for tabular output (not tracked by git)    
% Naming conventions: fields separated by "-" subfields by "_" 
    % Instances:    (with string $NAME and integer $SIZE)
        % $NAME_$SIZE-init.csv holds initial conditions; 
        % $NAME_$SIZE-trav.dat holds daily travelers matrix (air passengers +
        % commuters)
        %sample: a~NW~cty_75[-init.csv,-trav.dat], 75 counties in OR and WS
    % Output tables: $NAME_$SIZE-$CPG-[abs | frac]-out.csv
% 2021-03-04 v.0.1 up to reading IVs and travel matrix

%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is for problem instances (IVs + travel matrix)
instDir="data/by-tract";
oDirFig = "fig"; %write the output figures and tables here
if(~exist(oDirFig,"dir"))
    mkdir(oDirFig); %make sure it exists
end
oDirTab = "out"; %write the output figures and tables here
if(~exist(oDirTab,"dir"))
    mkdir(oDirTab); %make sure it exists
end

fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields
IV_suff="init.csv"; %all instances IVs end like this 
trav_suff="trav.dat"; %all instances' travel data end like this

%IV_Path=fullfile(instDir,instName+fnSep+IV_suff); %path to the IVs
pathIV= @(iname) fullfile(instDir,iname+fnSep+IV_suff); %path to the IVs
pathTrav = @(iname) fullfile(instDir,iname+fnSep+trav_suff); %path to the travel data, if any

%paths to output tables
pathOutTabAbs = @(iname) fullfile(outDir,join([iname,"abs"],fnSep));
pathOutTabFrac = @(iname) fullfile(outDir,join([iname,"frac"],fnSep));

%% Problem Instance (initial values, populations, and travel matrix)
%inst="a~NW~cty_75"; %by-county OR + WS, with flights & commute
inst="a~NW~ste_2"; %by-state OR + WS, with flights & commute

% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tIVs = readtable(pathIV(inst));

nodeNum = size(tIVs,1); %as many nodes as there are rows
N = table2array(tIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N etc.

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (Araw - diag(Araw) + diag(N));
%% Model Parameters
% these follow (El Ouardighi, Khmelnitsky, Sethi, 2020)

tFin = 180; %time is [0,tFin]