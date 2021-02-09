"""
Author: Yaroslav Salii, 2020+.

This is Oboe-Mangle, or maybe Arshyn, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS

Reading the FluTE data and transforming it to my input format.
Might separate my input format definitions later.

Oboe.jl v.0.1: "From scripts to proper code" edition
        v.0.2: switch to `module`
        v.0.3: got as far as reading FluTE tracts
        v.0.4: added basic by-county and by-state IV aggregation
        v.0.5: join (intersect) BTS with OpenFlights
        v.0.6: added designated AP assignment
'20-12-18   v.0.7: wrote pairs-to-matrix xform by hand (vs. `unstack`, which was unpredictable)
'20-12-30   v.0.8: added a beta node-to-node daily air passenger computation
'21-01-03   v.0.8.1: a half-baked Main(), look in bit bucket. Tested on 3K by-county!
'21-01-22   v.0.9: air travel and commute matrices are in, tested on 2K NW (by-tract)
'21-02-02   v.0.9.1: hotfix, removing deprecated `by`, `names!`, and the like (half-done)
'21-02-02   v.0.9.2: removed all the deprecated stuff, and it all works again
"""

#TODO: make debug defaults parameterized, via macros or otherwise

#module, not `include`, to prevent multiple inculdes (oh hi #ifndef)
#submodules may be includes though, when I get to carve it up
module Oboe


#reading ingress $name-tract.dat, $name-wf.dat$
#reading ingress; possibly, output too
using CSV
#transforming the data in tabular form;
using DataFrames
#ready-made `mean` function
using Statistics
#for `haversine` formula of distance between two points on a big circle
using Distances


#====BASE===FILENAMES==TYPES==DATA=STRUCTURES=====#

const callsign="This is Oboe v.0.9.2"
#println(callsign)

#=
1. all paths are set in view of running from epi-net-m/tools
2. data is meant to live in epi-net-m/data
=#

#all you need to know about input and output of instances
struct NamingSpec #all fields are String, don't say I didn't warn you
    # e.g. sep="-", sSep="_"; $name_$size-$suff
    sep::String #to the right of last $sep is filename suffix, to the left is the instance name
    sSep::String
    #read from ifDir, write to ofDir
    ifDir::String
 #   ifDirAir::String
    ofDir::String
    # what's after instance's name in its Initial Values file name
    fltInitSuff::String
    myInitSuff::String
end

#the naming conventions I am going to use
#as well as input and output directories
global const fn=NamingSpec("-","_"
    ,joinpath("..","data","by-tract","flute")
    ,joinpath("..","data","by-tract")
    ,"tracts.dat","init.csv")


#locating BTS and OpenFlights input files
const APdir= joinpath("..","data","by-tract","air")::String
#raw BTS data, with separate per-carrier flights
global const ifBTS=joinpath(APdir,"2019 BTS domestic.csv")::String
#raw OpenFlights AP data
global const ifAPs=joinpath(APdir,"Openflights airports.dat")::String
#TODO: consider wgetting from the original https://raw.githubusercontent.com/jpatokal/openflights/master/data/airports.dat


#===FluTE======READ=&=PROCESS===FluTE===TRACTS=============#
# sample usage: Oboe.lsTracts()[4] |> Oboe.rdFluteTract |> Oboe.aggBySte
# a *node* (tract) must have the fields :Pop,:LAT,:LNG

#show the FluTE's tract filenames found in ins.ifDir, default to fn
function lsTracts(ins::NamingSpec = fn)
    filter(s::String -> endswith(s,ins.fltInitSuff),readdir(ins.ifDir))
end

#= load a FluTE census tract info into a DataFrame,
setting the types and colnames 
Default to NW tracts (Oregon + Washington)
=#
function rdFluteTract(ifName::String=lsTracts()[findfirst(s -> startswith(s,"a~NW"),lsTracts())],
            ins::NamingSpec=fn)
    idf=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    rename!(idf, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
    idf.Name = map((s,z,w)-> join([s,z,w],"~"),idf.Ste,idf.Cty,idf.Tra)
    return idf
end

# read the FluTE's all-US-by-tract instance (usa-tracts.dat)
function rdWholeUS(ins::NamingSpec=fn)
    ifName = filter(s -> split(s,"-")[1]=="usa",lsTracts(ins) )[1]
    rdFluteTract(ifName,ins)
end

#--------AGGREGATE---TRACT-LIKE---DATA--------------#
#= by-state PRE-aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
TODO: put the by--end output into a variable and add post-processing (adding the ID)
=#
function aggBySte(idf=rdFluteTract()::DataFrame)
    gd = groupby(idf,[:Ste]) #group by U.S. State FIPS
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
end

#= by-county PRE-aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggByCty(idf=rdFluteTract()::DataFrame)
    gd = groupby(idf,[:Ste,:Cty]) #group by U.S. County FIPS, within the same State FIPS 
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
end

#=======WORKING===WITH===AIRPORTS====================#


#-------TIDY---OPENFLIGHTS---INPUT-------------------#

#move these consts to BASE?
const apAllColNames = [:ID,:Name,:City,:Country,:IATA_Code,:ICAO_Code,:LAT,:LNG,:Altitude,:Timezone,:Daylight_Savings,:TZ,:Type,:Source]
#reordered in the order of necessity; 3-letter IATA code as ID
const apRetainedColNames=[:IATA_Code,:LAT,:LNG,:Name,:City,:Country]

#=
read the OpenFlights.org's airports.dat, and retain only the columns
I feel I may use; put :IATA_Code,:LAT,:LNG first,
these are definitely useful
CAVEAT: some missing values in :City, some weird values in :Country
TODO: resolve this caveat
=#
function rdAPs(ifName=ifAPs::String)
    out=CSV.File(ifName,header=false) |> DataFrame
    rename!(out, apAllColNames)
    #retain only the useful columns
    select!(out,apRetainedColNames)
end


#============BTS===FLIGHT===DATA===PROCESSING============#

#---READ---BTS--FLIGHTS-------------------------#
#=read the BTS file, and retain only :1 Passengers, :5 ORIGIN, and :7 DEST
#set the columns to [:ORG,:DST,:PSG] for uniformity=#
function rdBTS(ifName=ifBTS::String)
    rawBTS = select(CSV.read(ifName,DataFrame), :1,:5,:7)
    if map( x -> floor(x)==x, rawBTS[:,1]) |> all #if :PSG are Integer
        rawBTS.PASSENGERS=convert(Array{Int64,1},rawBTS.PASSENGERS)
    else #∃ non-integer passengers-per-year entry
        error("CAVEAT: *non-integer* PASSENGERS in Flights input.")
    end
    rename!(rawBTS,[:PSG,:ORG,:DST])
    return select(rawBTS,[:ORG,:DST,:PSG]) #
end


#---BTS---TIDY--AGGREGATE---ETC---------------------------#

#=
sum the passengers on the flights with same (ORG,DST) pairs
sort lexicographically in (ORG,DST) order
final output is the air travel graph as a *list of edges*
=#
function grpBTS(idf=rdBTS()::DataFrame)
    gd = groupby(idf,[:ORG,:DST])
    out  = combine(gd,:PSG => sum, renamecols = false)
    sort!(out,[:ORG,:DST]) #sort by departure/arrival AP names
end


#input DF must have :ORG,:DST,:PSG cols
#CAVEAT: missings are set to 0.0
function mkFlightMx2(idf = grpBTS()::DataFrame; init_to = 0.0::Number)
    orgAPs= [ row.ORG for row ∈ eachrow(idf)]
    dstAPs = [ row.DST for row ∈ eachrow(idf)]
    allAPs = orgAPs ∪ dstAPs |> sort #make sure they are sorted by name (:IATA_Code in fact)
    dim = length(allAPs) #the matrix' will be [dim × dim]
    ixs = zip(allAPs, 1:dim) |> Dict #get index by name 
    xis = zip(1:dim,allAPs) |> Dict #get name by index
    #init the matrix with all init_to, defaulting to 0.0: screw the missings, I am nullifying them anyway
    outM = fill(init_to,(dim,dim))
    #fill the matrix with the *known* values
    for row ∈ eachrow(idf)
        outM[ixs[row.ORG],ixs[row.DST]] = row.PSG  
    end
    
    return (M= outM, ix=ixs,xi=xis,apCodes=DataFrame(IATA_Code=allAPs))
end

#=
in: apCodes [:IATA_Code], 
output: a `df` [:IATA_Code,:IN,:OUT,:TOUR],
with :IN=Σ_incoming PSG, :OUT=Σ_outgoing,:TOUR=Σ_(:ORG=:DST)
in :IN and :OUT sums, `missing` is non-absorbing and the diagonal is omitted
(opt) :TTL=:IN+:OUT,
=#
function mkAggFlows2(apCodes::DataFrame=mkFlightMx2().apCodes, M=mkFlightMx2().M::Matrix)
        out=DataFrame(IATA_Code=apCodes.IATA_Code
    #col-wise total sans the reflexive, `missing` if the arrivals are only reflexive
        ,IN=[ (sum(M[:,j]) == M[j,j]) ? missing : (sum(M[:,j]) - M[j,j]) for j ∈ 1:size(M)[2] ]
    #row-wise total sans the reflexive, `missing` if the departures are only reflexive
        ,OUT=[ (sum(M[i,:]) == M[i,i]) ? missing : (sum(M[i,:]) - M[i,i]) for i ∈ 1:size(M)[1] ]
    #no. reflexive travelers, or `missing`
        ,TOUR=[M[i,i] for i ∈ 1:size(M)[1]] )#just the reflexive travelers
    return out
end


#------BTS----CENSORING--------------------#

#TBD: enplanement is just OUT, so maybe not sum it with in.
#TBD: rename to dropSmallAPs
#only retain the APs with at least `p` annual enplanements (sum :IN and :OUT)
function censorAggFlows(p=2500::Number,idf=pickCleanAPs()::DataFrame)
    if (["IN","OUT"] ⊈ names(idf)) error("wrong DF") end
    filter(row -> row.IN+row.OUT ≥ p, eachrow(idf)) |> DataFrame
end

#throw out the APs with missing IN or OUT enplanements
scrubAPs(APs=mkAggFlows2()::DataFrame) = filter(row -> !ismissing(row.OUT) && !ismissing(row.IN)
,eachrow(APs)) |> DataFrame

#===========INTERCONNECT===================#

#=====BTS=<-->=OPENFLIGHTS=======#


#pick just the given APs (by IATA_Code) from all in dfAPinfo
#output columns as [:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country]
function pickAPs(myAPs=mkAggFlows2()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
    out= innerjoin(myAPs, dfAPinfo,on=:IATA_Code)
    select!(out,[:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country])
end

#= pick just the “clean”  APs (no `missing`) from all in dfAPinfo
CAVEAT: doesn't retain any APs that have :IN or:OUT missing or 0;
recall that 0 was identified with `missing` in mkAggFlows()=
Also remove missings from the :IN and :OUT
=#
function pickCleanAPs(myAPs=mkAggFlows2()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
    pickAPs(myAPs,dfAPinfo) |> scrubAPs
end

#----DESIGNATED---APs-----------#

#= Distance
Starting “haversine”/orthodromic/big-circle distance
Should get the almost true (ball vs. geoid)
when coupled with Earth's radius, or rather, AP's and loc's altitudes
CAVEAT: getting both altitudes is almost as bad 
as scrubbing “road distances” from the net
TODO: route through JuliaGeo/Geodesy
=#

#Earth's radius in km, IAU 2015
const R_Earth= 6371 
#x and y are Union{AbstractVector{T}, NTuple{2, T}} where T<:Real
myDist(x,y) = haversine(x,y,R_Earth)

#= Get Designated AP
for a node-row [:ID(:Ste,:Cty,:Tra),:Pop,:LAT,:LNG],
return the *nearest* AP's :IATA_Code and the distance :dst to it in km   

NB: worked slow-ish in Jupyter Notebook; might rework without sorting all this thing.
UPD: this slow-down was probably due to lack of caching with arguments by default
=#
function getDsgAP(node=aggBySte()[1,:]::DataFrameRow,APs=censorAggFlows()::DataFrame)
    alle=[(dst = myDist((node.LAT,node.LNG)
            ,(ap.LAT,ap.LNG))
            ,IATA_Code = ap.IATA_Code
            ,psg= ap.IN+ap.OUT)  for ap ∈ eachrow(APs)] 
    sort!(alle, by = first) # sort by distance to APs
    return (dsg=alle[1],all=alle)   
end
 
#= Assign a designated AP to each node in `nodes`, return a DF with info on both nodes and
their designated APs, [:ID,:Pop,:LAT,:LNG,:IATA_Code]; :IATA_Code for the designated AP
dbg columns: :dst distance to the dsg AP, :psg = :IN + :OUT passengers through dsg_AP
=#
function assignDsgAPs(nodes=aggBySte()::DataFrame,APs=censorAggFlows()::DataFrame)
    dsgAPs = map(n -> Oboe.getDsgAP(n,APs).dsg,eachrow(nodes)) |> DataFrame
    hcat(nodes,dsgAPs)
end

#----NODE-TO-NODE---PASSENGER--FLOWS----FROM--AP--DAILY--ENPLANEMENTS-----#
#CAVEAT: all nodes in this section must have a designated AP (:IATA_Code column)

#------------AUX-----------#
#say how many people fly through each AP; default to per-state aggregation
#`nodes` must have designated APs (:IATA_Code) and populations (:Pop)
function mkClusterPops(nodes=assignDsgAPs(aggBySte())::DataFrame)
    #group by the same designated AP: each group is this AP's __catchment area__
    gd = groupby(nodes,[:IATA_Code]) 
    #= the prime here is Σpopulation by designated AP, 
    but it's also a clusterization so here's `mean` :LAT and :LNG as dumb centroids =#
    out = combine(gd,:Pop => sum, :LAT => mean, :LNG => mean,renamecols = false)
end

#TODO: consider leaving it *inside* `assignDsgAPs`, it has to be called anyway
#make a `Dict` mapping each AP's ID (:IATA_Code) to the :Pop of its *catchment area*
function mkAP_pop_dict(nodes=assignDsgAPs(aggBySte())::DataFrame )
    cache=mkClusterPops(nodes)
    return Dict(cache.IATA_Code .=> cache.Pop)
end

#GLEaM-like aggregation: lump together all nodes in AP's catchment area (start with `census tracts`)
#difference with GLEaM: some airports may get lost, whereas GLEaM starts with APs (gotta re-check that!)
function aggByAPs()
    return assignDsgAPs(rdFluteTract()) |> mkClusterPops
end


#for a node `n`, the fraction of pop in `dsg_n`'s catchment area
#the node MUST have an :IATA_Code column (its *designated AP*)
function nodePsgShare(node = assignDsgAPs()[1,:]::DataFrameRow,dAP_pop=mkAP_pop_dict()::Dict)
    return node.Pop / dAP_pop[node.IATA_Code]
end

#apply the above to all nodes. 
#each node MUST have an :IATA_Code column (its *designated AP*)
function assignPsgShares(nodes=assignDsgAPs(aggBySte())::DataFrame,dAP_pop=mkAP_pop_dict()::Dict)
    psgShares = map(n -> Oboe.nodePsgShare(n,dAP_pop), eachrow(nodes)) #compute the shares
    out = hcat(nodes,DataFrame("shr" => psgShares)) #add them as a :shr column
end

#------AIR---PASSENGER---FLOW---MATRIX----------------------#
#NB! `nodes`:[:ID,:Pop,:IATA_Code,:shr]
#TODO: make usage of stopgap explicit, e.g. _(_;dbg=true) ... if dbg print("Debug info")
function mkPsgMx(ns=assignPsgShares()::DataFrame)
    retAPs = ns.IATA_Code |> unique #the APs that are designated for at least one `node`
    dim = nrow(ns) #final output matrix is [dim × dim]
    #retain only flights (tuples :ORG,:DST,PSG) for designated APs
    retFlows = filter( a -> a.ORG ∈ retAPs && a.DST ∈ retAPs, eachrow(grpBTS())) |> DataFrame
    #=find all the designated APs that didn't make it into `flowCns` 
    because they had no connections to other designated APs (“isolated dsg APs”) =#
    isolatedAPs = filter( apID -> apID ∉ retFlows.ORG && apID ∉ retFlows.DST, retAPs)
    #report if there were *any* isolated designated APs
    println("Found ", isolatedAPs |> length," ISOLATED designated APs: ",isolatedAPs)
    #=now make dummy, 0-passenger flights from each of them to the 1st :DST AP in retFlows
    to make sure `isolatedAPs` appear in `outM`, and add them to retFlows=#
    patchedFs = vcat(retFlows
            , [(ORG=apID,DST=retFlows[1,:].DST,PSG=0.0) for apID in isolatedAPs] |> DataFrame)
    sort!(patchedFs,[:ORG,:DST]) #restore the order (this screams for an object and a constructor!)
    #get the AP-to-AP flows for the designated APs, with IATA_Code <-> Index 
    aps = mkFlightMx2(patchedFs)  
    A = map(x -> x/365,aps.M) #make avg. daily ap-ap psg flows; perhaps just drop <1 values
    outM = fill(0.0, (dim,dim))
    #for each [from,to] pair, set 0.0 if dsg_APs match or weigh AP<->AP psg by nodes' pop shares
    for from ∈ 1:dim, to ∈ 1:dim
        if ns[from,:].IATA_Code ≠ ns[to,:].IATA_Code #no air travel unless dsg APs are different
            # shr_{from} × shr_{to} × psg_{dsg_from,dsg_to}
            outM[from,to] = ns.shr[from] *
                            ns.shr[to] *
                            A[ aps.ix[ns.IATA_Code[from]]
                                    ,aps.ix[ns.IATA_Code[to]] ]                                                        
        end
    end 
    return outM

end #end mkPsgMx()

#-------COMMUTER---FLOW--------------------#

# dim = length(allAPs) #the matrix' will be [dim × dim]
# ixs = zip(allAPs, 1:dim) |> Dict #get index by name 
# xis = zip(1:dim,allAPs) |> Dict #get name by index
# #init the matrix with all 0.0, screw the missings, I am nullifying them anyway
# outM = fill(0.0,(dim,dim))
# #fill the matrix with the *known* values
# [outM[ixs[row.ORG],ixs[row.DST]] = row.PSG  for row ∈ eachrow(idf)]

# return (M= outM, ix=ixs,xi=xis)


#read & tidy all "usa-wf-$fips.dat" for $fips in fipss; default to NW: Oregon + Washington
function rdTidyWfsByFIPS(fipss::Array{String,1}=["41","53"],ins::NamingSpec=fn)
    #generate FluTE filename by State's FIPS
    wf_by_FIPS(fips::String) = "usa-wf-$fips.dat"
    
    wfs = map(wf_by_FIPS, fipss)
    wfs2 = [CSV.File(ipath, 
        header = false,
        types =[String,String,String,String,String,String,Int64]) |> 
    DataFrame for ipath in map(f -> joinpath(ins.ifDir,f), wfs)]

    to_keep(r) = r[4] ∈ fipss && r[3]≠r[6]
    #TIDY: retain only the commute between states in `fipss`; this chucks the REFLEXIVE commute
    wfs3 = map(df -> filter(to_keep,eachrow(df) ) |> DataFrame,wfs2)
    wfs4 = (length(wfs3) > 1) ? reduce(vcat,wfs3) : wfs3 #add them all together
    #combine the State,County,Tract triples into single columns
    insertcols!(wfs4, 1, 
        ORG = map((s,z,w)-> join([s,z,w],"~"),wfs4.Column1, wfs4.Column2,wfs4.Column3))
    insertcols!(wfs4, 2, 
        DST = map((s,z,w)-> join([s,z,w],"~"),wfs4.Column4, wfs4.Column5,wfs4.Column6))
    rename!(wfs4, :Column7 => :CMT) #these are daily commuters between :ORG and :DST
    select!(wfs4,:ORG,:DST,:CMT) #chuck the unnecessary
    sort!(wfs4,[:ORG,:DST])

    return wfs4
end

#--------COMMUTER---FLOW---MATRIX--------------#


function mkCmtMx(ns::DataFrame,wfs::DataFrame)
    if "Name" ∉ names(ns) || !(["ORG","DST","CMT"] ⊆ names(wfs))
        error("Wrong DF! Terminating.")
    end
    dim = nrow(ns) 
    outM = fill(0.0,(dim,dim)) #commute matrix, node-to-node, ordered as `ns`
    ixs = zip(ns.Name, 1:dim) |> Dict #get index by name
    xis = zip(1:dim, ns.Name) |> Dict #get name by index
    for trip ∈ eachrow(wfs) #doggedly fill each trip into the matrix
        outM[ixs[trip.ORG],ixs[trip.DST]]=trip.CMT
    end #reflexive commute is removed beforehand

    return outM
end

#=========OUTPUT===PREP=========================#

#-----AUX----------------------#

#= Select an `id`-generating function
0. if it's a census tract, just write $Ste~$Cty~$Tra
1. if it's a county or a state, write as Integer to match `us-10m.json`
2. if it's a Voronoi cell of an AP, write the AP's IATA_Code
=#
function select_mkid(nms::Array{String})
    if ["Ste","Cty","Tra"] ⊆ nms #node's a census tract
         return r-> join([r.Ste,r.Cty,r.Tra],"~")
    elseif ["Ste","Cty"] ⊆ nms #node's a county, must accomodate `us-10m.json`
         return r-> parse(Int,r.Ste*"000") + parse(Int,r.Cty)
    elseif  ["Ste"] ⊆ nms #node's a state, must accomodate `us-10m.json`
         return r-> parse(Int,r.Ste)
     elseif ["IATA_Code"] ⊆ nms #node is Voronoi cell around an AP
         return r -> r.IATA_Code
     else
        error("Nodes header not recognized. Can't generate an :id without FIPS or IATA_Code.")
     end
 end

#= if the first node is sterile (I_i), seed it with 1 infected (Hi, I'm _idempotent_!)
input MUST have [:I_i,:S_i]=#
function infect!(ns::DataFrame)
    if ns.I_i[1] == 0
        ns[1,:S_i]-=1
        ns[1,:I_i]+=1
    end
    return ns
end
#-----INITIAL---VALUES----------------------#

#=
input MUST have [:IATA_Code,:Pop,:Name,:LAT,:LNG]
out spec: [:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG]
=#
function ns2iv_sterile(ns::DataFrame)
    out = copy(ns)
    rename!(out, :Pop => :N_i) 
    insertcols!(out,1,:id => map(select_mkid(names(out)), eachrow(out)))
    insertcols!(out,4, :S_i => out.N_i)
    insertcols!(out,5, :I_i => zeros(Int64,nrow(out)))
    insertcols!(out,6, :R_i => zeros(Int64,nrow(out)))
    select(out,[:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG])
end

#=========AUX---DBG==================#
#= This section has the functions I don't intend to move to production, which, however,
are useful enough to be pulled from disparate notebooks

Might leave them babbling out of stdout/stderr though
=#

#input: a square matrix of numbers; might add type-checking later
function talkDnsy(M)
    #no. ≠0 entries
    a=filter(x -> x> 0.0, M) |> length
    dim = size(M,1)
    b = dim*(dim-1) #max no. ≠0 connections
    dnsy = a/b
    println("We've got ",a," ≠0 connections")
    println("With at most ",b, " ≠0 entries, we get the density ",dnsy)
    #print("The density is ",filter(x -> x> 0.0, M) |> le )
    return(nnonzero = a, maxnonzero= b, density = dnsy, dim = size(M,1))
end

end #end module Oboe


#========BIT=====BUCKET===========#


#-----MAIN---CODE------------#
# #this piece is to be integrated through oboe-main.jl
# myshow = obj -> println(first(obj,5))
# Oboe.mkPsgMx(Oboe.assignPsgShares(Oboe.assignDsgAPs()))
# all AP-AP travel, as list o'pairs [:ORG,:DST,:PSG], :ORG and :DST are :IATA_Code
# pBTS = Oboe.grpBTS() 
# pBTS |> myshow
# #cache the smallest reasonable APs, [:IATA_Code,:LAT,:LNG,:IN.+:OUT ≥ 2500] 
# APs = Oboe.censorAggFlows() # 
# APs |> myshow
# #read the Northwest (OR,WA) census tracts
# nsRaw= Oboe.rdFluteTract("a~NW" * Oboe.fn.sep * Oboe.fn.fltInitSuff )
# nsRaw |> myshow
# #to each node, assign a designated AP from `APs` -> +[:IATA_Code]
# ns = Oboe.assignDsgAPs(nsRaw,APs) 
# #first(ns,5) |> print
# ns |> myshow
# #now find the :Pop of each APs' catchment area, and chuck that into a `Dict`
# d = Oboe.mkAP_pop_dict(ns) #406 APs end up `designated` for aggByCty, not bad
# #go on, compute the nodes' passenger shares (add the :shr col), with `d` in mind
# ns2 = Oboe.assignPsgShares(ns,d)
# ns2 |> myshow
# #finally, compute NODE-NODE daily air passengers
# nnPsg=Oboe.mkPsgMx(ns2)


# fipsNW=["41","53"]
# iwfs=Oboe.rdTidyWfsByFIPS(fipsNW)
# nnCmtMx = Oboe.mkCmtMx(ns2,iwfs)