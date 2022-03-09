%% Compute the Dynamics of Net-SIR with Fixed Control.
% Better IVs by running null-control pandemic

%Run the Net-SIR model with null control 
%to determine when does infection level reaches 1 in 10^5 everywhere
%and use it as initial conditions for the optimal control problem

%Author: Yaroslav Salii, 2022+
% Directories: 
    %./data/by-tract is for problem instances (IVs + travel matrix)
    %./fig is for the figures (not tracked by git)
    %./out is for tabular output (not tracked by git)    
% Naming conventions: fields separated by "-" subfields by "_" 
    % Instances:    (with string $NAME and integer $SIZE)
        % $NAME_$SIZE-init.csv holds initial conditions; 
        % $NAME_$SIZE-trav.dat holds daily travelers matrix (air passengers +
        % commuters)
        %sample: NWcty_75[-init.csv,-trav.dat], 75 counties in OR and WS
    % Output tables: as input, to /data/inst-100K
% 2022-02-23 v.0.0 chiseling away the unnecessary
% 2022-02-28 v.0.1 first and last past threshold rough-in
% 2022-03-01 v.0.1.1 ignore isolated, sterile nodes in thresholds (set 0)
%            v.0.1.2 abs + rel thresholding ( > nrm=10^-5 && > 1 person)
% 2022-03-02 v.0.2 new IVs from start threshold at selected node
% 2022-03-09 v.0.3 threshold IVs start at the min. time + steady state out
%% TODO
%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is for problem instances (IVs + travel matrix)
instDir="../data/by-tract";
ofigDir = "../fig"; %write the output figures and tables here
if(~exist(ofigDir,"dir"))
    mkdir(ofigDir); %make sure it exists
end
otabDir = "../data/inst-100K-m1"; %write the output tables here
if(~exist(otabDir,"dir"))
    mkdir(otabDir); %make sure it exists
end

fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields


%YS: using regexp/extract as opposed to split by fnSep="-" because
%it's too complicated to extract the first element of split's return array
%in an anonymous function

%instance name as regexp pattern; 
rpInst = regexpPattern('[A-Z]+[0-9]*(tra|cty|ap|ste)_[0-9]+');

IV_suff="init.csv"; %all instances IVs end like this 
trav_suff="trav.dat"; %all instances' travel data end like this

%INPUT PATHS & FILENAMES
pathIV= @(iname) fullfile(instDir,iname+fnSep+IV_suff); %path to the IVs
pathTrav = @(iname) fullfile(instDir,iname+fnSep+trav_suff); %path to the travel data, if any

%OUTPUT PATHS & FILENAMES
pathNewIV = @(iname) fullfile(otabDir,iname + fnSep + IV_suff);
%LOG FILE, a CSV with a line for each iteration
pathotlog= @(iname) fullfile(otabDir,iname+"-log.csv"); 

%% Model Parameters
% these follow (El Ouardighi, Khmelnitsky, Sethi, 2020)

%INFECTION and RECOVERY RATES
%(Rader, Scarpino, Nande, &al., 2020): alpha = 0.15, gamma = 0.1
%baseline R_0=2.74; beta/gamma 63%
%set to harmonic mean of \alpha and \beta from idem
beta = 0.1196; % 1/(1/ 0.2977 + 1/0.2); %infection rate; (idem)
%1/beta = 8.36 mean time to be (symptomatic) infected
%beta = 1/2.5; %compat with older, R_0 = 3.3
gamma  = 0.0437; % 1/ (beta/R_0); R_0 = beta / gamma (idem)
%1/gamma = 22.904 mean time to recovery
%gamma = 1/8.3; %compat with older

%THRESHOLDS etc.

nrm = 10^-5; % target 1 in 100K infected
%


%RUNNING COSTS      l for lockdown (control)
c = 0.00200; %running cost of infections; mean(c_1=100, c_2=300, c_3=200)
l = 0.000450; %running cost of control; (q3=450, securing social interactions)
%TERMINAL COST
k = 0.02000; %terminal cost of infections; mean(k_1=1K, k_2=3K, k_3=1K)
%FATIGUE RATES      (PI/PD-SF) 
r1 = 0.002; %infection rates fatigue rate
%r1 = 0; %infection rates fatigue rate
r2 = 0.002; %lockdown control fatigue rate
%r2 = 0;

%%TIME
T = 1000; %time is [0,T], in days; default to 500 to hopefully reach endemic
%tArr = 0:T that's the de facto usage

%CONTROL bounds (for each component, at each time)
umin = 0; umax = 1; % u_i \in [0,1] \forall i \in nodes

%% Logs Setup
% a table to log dynamics computation times etc.
tablogNames = ["Inst","t_hms","t_ms",...
    "pop","cZR","RR","ZZ"];

tablogTypes = ["string","duration","double",...
    "double","double","double","double"];
    
tablogSize=[1 size(tablogNames,2)];
tablog = table('Size',tablogSize,'VariableType',tablogTypes,'VariableNames',tablogNames);

% a table to log the node chosen to fix the start of the epidemic
% see its signature in "Determine the start day" section
tabsummary = table;

%% Set Problem Instance Name

% read all instances IV files names' into array-of-structs
sAllIVs = dir(fullfile(instDir,'*-init.csv')); %array o' structs
cellAllIVs = extractfield(sAllIVs,'name'); %all IV file names (cell array)
cellAllInst = cellfun(@(s) extract(s,rpInst),cellAllIVs); %all instance names

%inst="NWtra_2072"; %by-tract OR + WS, with flights & commute
inst="NWcty_75"; %by-county OR + WS, with flights & commute
%inst="NWap_23"; %by-airport OR + WS, with filghts & commute
%inst="NWste_2"; %by-state OR + WS, with flights & commute

%inst="WCTtra_9110"; %WCT is CA + OR + WS
%inst="WCTcty_133";
%inst="WCTap_52";
%inst="WCTste_3";

%inst="ALLcty_3109";
%inst="ALLap_417";
%inst="ALLste_49";

%inst = "NW1tra_2834";
%inst = "NW1cty_136";
%inst = "NW1ap_37";
%inst = "NW1ste_4";

%inst = "NW2ste_8";
%inst = "NW2ap_74";
%inst = "NW2cty_259";
%inst = "NW2tra_4830";

%inst = "CAste_1"; %implementation doesn't support 1-node instances :(
%inst = "CAap_32";
%inst = "CActy_58";
%inst = "CAtra_7038";


for inst = cellAllInst(15:15) %normally, no.16 is "NWste_2"
    disp(pathIV(inst));
%end

%% Read Problem Instance (initial values, populations, and travel matrix)
% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tabIVs = readtable(pathIV(inst));

n = size(tabIVs,1); %as many nodes as there are rows
N = table2array(tabIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N (--> A)

NN = N / max(N); %normalized to max pop---for J, etc.
iNN = arrayfun(@(x) 1/x,NN); %inverse norm-pops, for Hadamard division (--> lax,u,J)
%NN = ones(n,1); iNN = ones(n,1); %FAT CAVEAT: "egalitarian" pricing

%STATE: initial conditions
s0 = table2array(tabIVs(:,4)) .* iN; %susceptibles at t=0, frac
z0 = table2array(tabIVs(:,5)) .* iN; %infecteds at t=0, frac
x0 = [s0;z0]; %overall state is [s;z]

%COSTATE: terminal values from transversality conditions (column vectors!)
%unused but I may find a use for these later
lasT = zeros(n,1); %\lambda_s(T) = 0_n
lazT = exp(r1*T)*k*NN; %\lambda_z(T) = e^{r_1T}k\tilde{N}
laxT = [lasT;lazT]; %costate is [\lambda_s; \lambda_z] 

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
A = Araw'; %easier to add one transpose here than in every equation
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (A - diag(A) + diag(N));
%A = eye(n); %dumb debug: isolated nodes
%A = ones(n,n); %dumber debug: homogeneous mixing

%% Population Stats
 %Nq20 = quantile(N,20);  Nq4 = quantile(N,4); 
 Nq10 = quantile(N,10);


%% Run Initialization

%STATE: 0-init [2n \times |tArr|] + initial conditions
x = zeros(n*2,T+1); x(:,1) = x0;
%COSTATE: 0-init [2n \times |tArr|] + terminal conditions
lax = zeros(n*2,T+1); 
lax(:,end) = laxT; 

%CONTROL
u = zeros(n,T+1); %constant zero control
ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines

%MISC
% delta = 1E-3; %min. relative error for norms of u,x,lax in stopping conditions
% ct = 0; %set the loop counter
opts = odeset('InitialStep',1,'AbsTol',1e-9,'RelTol',1e-6); %ensure the 1st step is at most 1 day long

%Run Flags
flg.computeCostate = 0;
flg.computeObjval = 0;
%% Compute the Dynamics
tNumdynStart = tic;    
    
    %STATE
    ftx1 = @(t,x) futxp(ppval(ppu,t),t,x,beta,gamma,A);
   % disp('Run ode45 on IVP for state x---forwards from 0 to T')
    x_sln = ode45(ftx1, [0 T], x(:,1),opts);
    x = deval(x_sln, 0:T );
    tStateDynDone = toc(tNumdynStart); tHMS = seconds(tStateDynDone);
    
    %COSTATE
    if (flg.computeCostate)
        gtx1 = @(t,lax) guxtlp(ppval(ppu,t), deval(x_sln, t) ...
            , t, lax, beta, gamma, A, r1, c, NN);
    %   disp('Run ode45 on IVP for costate lax---backwards from T to 0');
        lax_sln = ode45(gtx1, [T 0], lax(:,end), opts);
        lax = deval(lax_sln, 0:T);
    end
    
    %OBJECTIVE FUNCTION
    if (flg.computeObjval)
        Lt1 = @(tArr) Ltxu(ppval(ppu,tArr),deval(x_sln,tArr),tArr,r1,r2,c,l,NN); %running cost
     %   disp('Compute the objective function J(u,x,T)');
        PsiT1 = PsiT(x(:,end),T,r1,NN,k); %terminal cost
        J = PsiT1 + integral(Lt1,0,T); %the objective function
    else
        J = 0;
    end
    %STATISTICS
    
    tStatStart = tic;
    %recovered at day T
    RR_T = sum(N - (x(1:n,end) + x(n+1:end,end) ) .* N); 
    %infected at day T + recovered at day T: eff. all but susceptible at T
    cZR_T = sum( (ones(n,1) - x(1:n,end)) .* N);
    % infected at day T
    ZZ_T = sum(x(n+1:end,end) .* N);
    
    %cumulative infected, including recovered, daily
    %repmat(N, [1,numel([0:T])] ) - x(1:n,:).* N
    
    tStatDone = toc(tStatStart);
    
    %LOG
    currLog = {inst,tHMS,tStateDynDone,...
    sum(N),round(cZR_T,4),round(RR_T,4),round(ZZ_T,4)};
    
    if (ismissing(tablog.Inst(1))) %detect if it's the first iteration
        tablog(1,:) = currLog; %overwrite the dummy entry
    else
        tablog = [currLog; tablog]; %add the new entry on top of the log
    end
    disp(tablog(1,:)); %show the current log
    
%tNumdynEnd = toc(tNumdynStart);

%% Determine the start, peak, and end times for each node

%slice the state into (s,z,r) compartments
s = x(1:n,:); z = x(n+1:end,:); r = (1 - s - z); %[s z r] for output
S = s .* N; Z = z .* N; R = r .* N; cZR = Z + R;

%try thresholding on absolute people > 1
fndfirstabovethr2 = @(arrRel,node,thr,arrAbs) ...
    find(arrRel(node,:) > thr & arrAbs(node,:) > 1 ,1);
ffanrm2 = @(node) maybevalue(fndfirstabovethr2(z,node,nrm,Z),double(0));
fndlastabovethr2 = @(arrRel,node,thr,arrAbs) ...
    find(arrRel(node,:) > thr & arrAbs(node,:) > 1 ,1,'last');
flanrm2 = @(node) maybevalue(fndlastabovethr2(z,node,nrm,Z),double(0));

%search in all nodes (thresholds are in the functions)
nodeSet = 1:size(Z,1);

%find days when infected first go above and first go below the threshold
arrFirstAbove2 = arrayfun(ffanrm2, nodeSet);
arrLastAbove2 = arrayfun(flanrm2, nodeSet);

%find the peak day (nb! peakVal in non-normalized, absolute persons)
[arrPeakVal, arrPeakDay] = max(Z');
%disp(arrFirstAbove');

%some descriptive-like stats on the outbreak
tabdesc = table;

%consider adding tIVs.id, .Name
tabdesc.id = tabIVs.id;
tabdesc.pop = N;
tabdesc.peak = arrPeakVal';
tabdesc.peakFrac = arrPeakVal' ./ N;
tabdesc.start  = arrFirstAbove2';
tabdesc.peakDay = arrPeakDay'; 
tabdesc.end = arrLastAbove2';
tabdesc.Name = tabIVs.Name;

disp(tabdesc);


%% Determine the thresholds and peak in the *overall* model (as if it was one huge node)
ovrPop = sum(N); Z_ovr = sum(Z,1); z_ovr = Z_ovr ./ ovrPop;
% get the overall thresholds (when summed to single-node model)
ovrFirstAbove = find(z_ovr > nrm & Z_ovr > 1,1,'first');
ovrLastAbove = find(z_ovr > nrm & Z_ovr > 1,1,'last');
% get the overall peak
[ovrPeakVal, ovrPeakDay] = max(Z_ovr); 
ovrPeakFrac = z_ovr(ovrPeakDay);
cZR_ovr = sum(cZR(:,end));
ovrAffectedFrac = cZR_ovr / ovrPop;

tabovrNames = ["ovrStart","ovrPeakDay","ovrEnd","ovrPeak","ovrPeakFrac","ovrPop","ovr_cZR","ovr_cZR_%",];
tabovrTypes = ["uint32","uint32","uint32","double","double","uint32","uint32","double"];
tabovrSize = [1 size(tabovrNames,2)];
tabovr = table('Size',tabovrSize,'VariableType',tabovrTypes,'VariableNames',tabovrNames);
tabovr(1,:) = {ovrFirstAbove,ovrPeakDay,ovrLastAbove,round(ovrPeakVal,4),round(ovrPeakFrac,4),ovrPop,ceil(cZR_ovr),round(ovrAffectedFrac,4)*100};
% 
%% Determine the start day
% pick the *start* day of a *meaningful* node 
% = above Nq10(1) if min(N) < 4000 (arbitrary)
% = start day is above 40
% = no need to control for the *seed* node when t > 40

% set minimal node pop above first decile if minimal population is < 4000
popThreshold = 0;
if (min(N) < 4000)
    popThreshold = Nq10(1); %first decile
end

% ensure early starts go first
tabdescByStart = sortrows(tabdesc,{'start'},{'ascend'});
%select the benchmark node
ixbNode = find(tabdescByStart.start > 40 & tabdescByStart.pop > popThreshold,1);
% ixbNode is a node's index; now copy its record and add the instance name
bNodeEntry = tabdescByStart(ixbNode,:);
rowsummary = addvars(bNodeEntry,inst,'NewVariableNames',{'Inst'},'Before','id');
%add the overall info as well
rowsummary = [rowsummary,tabovr(1,:)]; %
disp("Fixing epidemic start at start day of the following node:");
disp(rowsummary);%disp(tabdesc(ixbNode,:));

%add to the summary if it wasn't empty
if(isempty(tabsummary))
    tabsummary = rowsummary;
else
    tabsummary = [rowsummary; tabsummary];
end

%% Make and write new IVs based on the state at tabsummary.start(1)

tabnewIVs = tabIVs;
tabnewIVs.S_i = S(:,rowsummary.start);
tabnewIVs.I_i = Z(:,rowsummary.start);
tabnewIVs.R_i = R(:,rowsummary.start);

% write the new initial values
%writetable(tabnewIVs,pathNewIV(inst));

end %end for the big instance-wise loop

%% Done inst-wise processing, work with overall info


% write the summary
%writetable(tabsummary,fullfile(otabDir,"compdyn-summary.csv"));

% write the solution output tables
% writetable(otabs,pathotabs(inst));
% writetable(tablog,pathotlog(inst));

%% AUXILIARY FUNCTIONS

%Return 0 when isempty or the value if it's not
function maybevalue = maybevalue(val,emptysigil)
    if(~isempty(val))
        maybevalue = val;
        return
    else
        maybevalue = emptysigil;
        return 
    end
end

%Hi, I'm the terminal cost in the objective J and I'm too small for a separate file
function PsiT = PsiT(xT,T,r1,NN,k)
n = size(xT,1) / 2; zT = xT(n+1:end,:);
PsiT = exp(r1*T) * k * dot(zT,NN);
end


%% BIT BUCKET

%normalize to 1-node model: sum all absolutes, and divide by pop
% A1 = @(c) sum( c .* N ); 
% a1 = @(c) A1(c) / sum(N);
% make the absolutes
% S = s .* N; Z = z.*N; R = r.*N; %multiply to get the absolute state
% S01 = sum(sNull .* N); Z01 = sum(zNull .* N); R0 = sum(rNull .* N);
% s11 = a1(s); z11 = a1(z); r11 = a1(r); 
% s01 = a1(sNull); z01 = a1(zNull); r01 = a1(rNull);

%column names for frac and abs tables (will be set to uppercase for abs)
% cns = [arrayfun( @(n) 's'+string(n),0:T) ...
%      arrayfun( @(n) 'z'+string(n),0:T) ...
%      arrayfun( @(n) 'r'+string(n),0:T)];
% 
% %cns_u = [ arrayfun( @(n) 'u'+string(n),0:T) ]; %col names for controls
% %cns_frac_u = horzcat(cns,cns_u); %combined frac with control
% 
% %col. names for per-node average infected-null, infected-opt, control effort
% cns_avgc = ["zNull_avg","z_avg","u_avg"];
%  
% %construct the solution output tables
% otabs = horzcat(tIVs,array2table([s z r] .* N,'VariableNames',upper(cns)));
% 
% %construct the average output table
% otabavgc = array2table(avgOut, 'VariableNames',cns_avgc);
% %note that the log table is already constructed
