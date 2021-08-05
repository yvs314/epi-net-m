"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

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
'21-02-14   v.0.9.3: add node partitioning, streamline mkFlightMx2() and mkPsgMx2() with aux
'21-02-14   v.0.9.4: add partition-to-partition flight matrix
'21-02-18   v.0.9.5: add partition-to-partition commute matrix
'21-03-31   v.0.9.6: min,median,max for nonzero elements in talkDnsy
            v.0.A: add by-AP agg with recomputation bypass in mkPsgMx()
'21-07-14   v.0.A.1: add checkIfFilesExist and censorFluteTractByFIPS
'21-07-19   v.1.0: fixed issue where unmatched tract IDs led to a crash
"""

#TODO: make debug defaults parameterized, via macros or otherwise

#module, not `include`, to prevent multiple inculdes (oh hi #ifndef)
#submodules may be includes though, when I get to carve it up
module Oboe



#CSV: IO FluTE & my, DelimitedFiles: writing the matrices (.dat)
using CSV,DelimitedFiles
#transforming the data in tabular form;
using DataFrames
#ready-made `mean` function
using Statistics
#for `haversine` formula of distance between two points on a big circle
using Distances
#for `diag` in `talkDnsy`
#using LinearAlgebra


#====BASE===FILENAMES==TYPES==DATA=STRUCTURES=====#

const callsign="This is Oboe v.1.0"
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

function checkIfFilesExist(iname::String; ofdir=fn.ofDir)
    dir = readdir(ofdir)
    any(filename -> startswith(filename, iname), dir)
end

function writeMe(iname::String,ivs::DataFrame,A::Array{Float64,2};ofdir=fn.ofDir)
    ofivs= join([iname,nrow(ivs)],"_") * "-init.csv"
    oftrv= join([iname,nrow(ivs)],"_") * "-trav.dat"
    CSV.write(joinpath(ofdir,ofivs),ivs)
    writedlm(joinpath(ofdir,oftrv),A)
end


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

# given a DataFrame of tracts, remove those that aren't from the given states
function censorFluteTractByFIPS(tracts::DataFrame, fips::Vector{String})
    filter(row -> row.Ste ∈ fips, tracts)
end



#--------AGGREGATE---TRACT-LIKE---DATA--------------#

 #NB! :Name is a String
 function select_mkName(nms::Array{String})
    if ["Ste","Cty","Tra"] ⊆ nms #node's a census tract
         return r-> join([r.Ste,r.Cty,r.Tra],"~")
    elseif ["Ste","Cty"] ⊆ nms #node's a county
         return r-> join([r.Ste,r.Cty],"~")
    elseif  ["Ste"] ⊆ nms #node's a state
         return r-> string(r.Ste) #ensure it's output as a string
     elseif ["IATA_Code"] ⊆ nms #node is Voronoi cell around an AP
         return r -> r.IATA_Code
     else
        error("Nodes header not recognized. Can't generate a :Name without FIPS or IATA_Code.")
     end
 end

#= by-state PRE-aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
TODO: put the by--end output into a variable and add post-processing (adding the ID)
=#
function aggBySte(idf=rdFluteTract()::DataFrame;make_names=true)
    gd = groupby(idf,[:Ste]) #group by U.S. State FIPS
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
end

#= by-county PRE-aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggByCty(idf=rdFluteTract()::DataFrame;make_names=true)
    gd = groupby(idf,[:Ste,:Cty]) #group by U.S. County FIPS, within the same State FIPS 
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
end

#= by-AP (:IATA_Code) PRE-aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat WITH dsg AP, a la [:Name,:IATA_Code,:Pop,:LAT,:LNG]
=#
function aggByAP(idf::DataFrame;make_names=true)
    gd = groupby(idf,[:IATA_Code]) #group by the tract's dsg AP :IATA_Code 
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
    sort!(out,:IATA_Code) #ensure it's sorted by the name/IATA_Code
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


#=
Transform a *list of arcs* into *adjacency matrix*
input: DF must have [:ORG,:DST,:PSG] cols
=#
function mkFlightMx2(fs=grpBTS()::DataFrame; daily=false::Bool,babble=false::Bool, init_to=0.0::Number)
    #make a sorted list of ALL the APs in `fs`
    allAPs = sort( (fs.ORG |> unique) ∪ (fs.DST |> unique))
    #delegate to the method with EXPLICIT ground set (`iretAPs`); with retaintours = true
    return mkFlightMx2(fs, allAPs; daily=daily,babble=babble,init_to=init_to,retaintours=true)
end

#takes *explicit* list of vertices `iretAPs`
function mkFlightMx2(fs::DataFrame, iretAPs::Array{String}; 
    daily=true::Bool # divide by 365, the input was annualized
    ,retaintours = false::Bool #flights with ORG==DST not retained by default
    ,babble=true::Bool, init_to=0.0::Number)
    
    retAPs = sort(iretAPs) #ensure AP codes are sorted, for indexing porposes
    pret = a -> a.ORG ∈ retAPs && a.DST ∈ retAPs
    if !retaintours
        pret = a -> a.ORG ≠ a.DST &&   a.ORG ∈ retAPs && a.DST ∈ retAPs
    end
    #pretain(a)  ? (retaintours) : ()
    #retain only flights (tuples :ORG,:DST,PSG) for APs ∈ retAPs; 
    retFlows = filter(pret, eachrow(fs)) |> DataFrame
    
    if babble #report if there were *any* isolated APs
        isolatedAPs = filter( a -> a ∉ retFlows.ORG && a ∉ retFlows.DST, retAPs)
        println("Found ", isolatedAPs |> length," isolated APs: ",isolatedAPs) 
    end
    
    dim = length(retAPs) 
    ix_ = zip(retAPs, 1:dim) |> Dict #get index by name 
    xi_ = zip(1:dim,retAPs) |> Dict #get name by index 
    M_ = fill(init_to,(dim,dim)) #default is 0.0; screw `missing`
    #fill the matrix with the *known* values
    for row ∈ eachrow(retFlows)
        M_[ix_[row.ORG],ix_[row.DST]] = row.PSG  
    end


    outM = daily ? map(x -> x/365,M_) : M_ #annual-to-daily, if requested
    return (M=outM, ix = ix_, xi = xi_,apCodes = retAPs)
end

#=
in: apCodes [:IATA_Code], 
output: a `df` [:IATA_Code,:IN,:OUT,:TOUR],
with :IN=Σ_incoming PSG, :OUT=Σ_outgoing,:TOUR=Σ_(:ORG=:DST)
in :IN and :OUT sums, `missing` is non-absorbing and the diagonal is omitted
(opt) :TTL=:IN+:OUT,
=#
function mkAggFlows2(apCodes::Array{String}=mkFlightMx2().apCodes, M=mkFlightMx2().M::Matrix)
        out=DataFrame(IATA_Code=apCodes
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
#NB! `ns`:[:ID,:Pop,:IATA_Code,:shr]
#NB! `outindices` must be at least the same length as `nrow(ns)`
#NB! All indices in `outindices` must not exceed `dim`
#This function calculates the air travel flow between each node pair in `ns`
#and accumulate it in a dim × dim adjacency matrix M, at a position determined
#by the `outindices` entries whose indices corresponds to the sequential indices
#of the nodes in the pair. So M[i1, i2] is the total passenger flow 
#from all nodes that have the `outindex` i1
#to all nodes that have the `outindex` i2.
function mkPsgMx(ns::DataFrame,
                 outindices::AbstractArray{Int},
                 dim::Int,
                 accumulate::Bool=true)
    insertcols!(ns, :outindex => outindices)
    #group the nodes into subDFs by their designated airport
    groupedbyAP = groupby(ns, :IATA_Code)
    numAPs = length(groupedbyAP)
    #from this, get list of IATA codes of unique airports that have nodes assigned
    retAPs = ns.IATA_Code |> unique
    #get the daily AP-to-AP flows for the designated APs, with IATA_Code as index
    fmx = mkFlightMx2(grpBTS(),retAPs; daily=true)
    #initialize output matrix
    outM = fill(0.0, (dim,dim))
    #for each combination of airports (AP1, AP2)...
    #(this weird iteration is to avoid repeating pairs)
    for i1 ∈ 1:numAPs, i2 ∈ (i1 + 1):numAPs
        # get IATA codes of AP1 and AP2
        ap1 = retAPs[i1]
        ap2 = retAPs[i2]
        # get AP1->AP2 and AP2->AP1 passenger numbers
        flow1to2 = fmx.M[fmx.ix[ap1], fmx.ix[ap2]]
        flow2to1 = fmx.M[fmx.ix[ap2], fmx.ix[ap1]]
        # skip over airport pairs that have no travel between them
        if flow1to2 == 0.0 && flow2to1 == 0.0
            continue
        end
        #for each pair of nodes (n1 ∈ AP1, n2 ∈ AP2)...
        grp1 = groupedbyAP[i1]
        grp2 = groupedbyAP[i2]

        indices1 = collect(grp1.outindex) :: Vector{Int}
        indices2 = collect(grp2.outindex) :: Vector{Int}
        shares1  = collect(grp1.shr)      :: Vector{Float64}
        shares2 =  collect(grp2.shr)      :: Vector{Float64}

        for n1 ∈ 1:nrow(grp1), n2 ∈ 1:nrow(grp2)
            #calculate psg flow n1->n2 and n2->n1 and store the result
            outindex1 = indices1[n1]
            outindex2 = indices2[n2]
            #do not count travel within the same subdivision
            if outindex1 == outindex2
                continue
            end
            shr1 = shares1[n1]
            shr2 = shares2[n2]

            if accumulate
                outM[outindex1, outindex2] += psg(shr1, shr2, flow1to2)
                outM[outindex2, outindex1] += psg(shr2, shr1, flow2to1)
            else
                outM[outindex1, outindex2] = psg(shr1, shr2, flow1to2)
                outM[outindex2, outindex1] = psg(shr2, shr1, flow2to1)
            end
        end
    end
    select!(ns, Not(:outindex)) #drop the outindex column from ns just in case
    return outM
end

function psg(fromshr, toshr, totalflow)
    fromshr * toshr * totalflow
end

#NB! `ns`:[:ID,:Pop,:IATA_Code,:shr]
#Returns an adjacency matrix representing the air passenger flow
#between all nodes in `ns`.
function mkPsgMx(ns::DataFrame)
    dim = nrow(ns)
    mkPsgMx(ns, 1:dim, dim, false)
end

#=
NB! now a method for partition-to-partition flights
`ns` MUST have [:IATA_Code,:shr]; `pns` MUST have [:Name]; `prt` :Name => [ns_row_indices]
=#
function mkPsgMx(ns::DataFrame,pns::DataFrame,prt::Dict;force_recompute=false)
    #first, assuming force_recompute is disabled, try to detect by-AP aggregation...
    retAPs = ns.IATA_Code |> unique |> sort
    if !force_recompute && "IATA_Code" ∈ names(pns) &&
        pns.IATA_Code == retAPs
        #...and if it is, just return the by-AP travel matrix
        println("By-AP aggregation engaged. Using AP-to-AP travel matrix directly.")
        return mkFlightMx2(grpBTS(), retAPs; daily=true).M
    end
    #build a dict mapping the name of a partition to its index in pns
    partindices = ((name, findfirst(n -> n == name, pns.Name)[1]) for name in keys(prt)) |> Dict{String, Int}
    #build an array mapping each node in ns to the index of its partition in pns
    outindices = zeros(Int, nrow(ns))
    for partname ∈ keys(prt), nodeindex ∈ prt[partname]
        outindices[nodeindex] = partindices[partname]
    end

    mkPsgMx(ns, outindices, length(prt))
end

#---AUX::---PARTITION---AND---REVERSE------------------------#
function partByCty(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.Ste[ri]==split(x,"~")[1] && 
                        ns.Cty[ri]==split(x,"~")[2],
                        1:nrow(ns)) for x ∈ pns.Name)
end

#consider a unified approach to partBySte and partByAP;
function partBySte(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.Ste[ri]==x, 1:nrow(ns)) for x ∈ pns.Ste |> unique)
end

function partByAP(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.IATA_Code[ri]==x, 1:nrow(ns)) for x ∈ pns.IATA_Code |> unique)
end


#part:: Name1 => [ns_row_indices],
#out:: ns_row_index => Name1 
function revexplPart(dict::Dict)
    ( v => k  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

#like above but also restore the Names of entries
# ns MUST have [:Name]
function revexplPart(dict::Dict, ns::DataFrame)
    ( ns.Name[v] => k  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

#like above but also map the keys' names to `indices`
# ns MUST have [:Name]; ix:: Name => partition_index
function revexplPart(dict::Dict, ns::DataFrame,ix::Dict)
    ( ns.Name[v] => ix[k]  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

#======COMMUTER=====FLOW==================#

#--AUX:--READ--CMT------#
#show the wf names
function lsWf(ins::NamingSpec = fn)
    filter(s::String -> startswith(s,"usa-wf-"),readdir(ins.ifDir))
end

ls_fipsAll = map(s -> split(split(s,".")[1],"-")[3] |> string,lsWf()) 


#read & tidy all "usa-wf-$fips.dat" for $fips in fipss; default to NW: Oregon + Washington
#if the tracts argument is provided, commutes referencing tracts that aren't in it
#will be discarded
function rdTidyWfsByFIPS(fipss::Array{String,1}=["41","53"], tracts::DataFrame=DataFrame(),
     ins::NamingSpec=fn)
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
    wfs4 = (length(wfs3) > 1) ? reduce(vcat,wfs3) : wfs3[1] #add them all together
    #combine the State,County,Tract triples into single columns
    insertcols!(wfs4, 1, 
        :ORG => map((s,z,w)-> join([s,z,w],"~"),wfs4.Column1, wfs4.Column2,wfs4.Column3))
    insertcols!(wfs4, 2, 
        :DST => map((s,z,w)-> join([s,z,w],"~"),wfs4.Column4, wfs4.Column5,wfs4.Column6))
    rename!(wfs4, :Column7 => :CMT) #these are daily commuters between :ORG and :DST
    select!(wfs4,:ORG,:DST,:CMT) #chuck the unnecessary
    #Chuck commute pairs that have tract IDs not present in tracts
    wfs5 = isempty(tracts) ? # check if tracts has been passed
        wfs4 :
        scrubWfs(tracts, wfs4) |> DataFrame
    sort!(wfs5,[:ORG,:DST])

    return wfs5
end

#=return wfs excluding the commutes for tracts not present in `ns`
ns MUST have [:Name]; wfs MUST have [:ORG,:DST,:CMT]
=#
function scrubWfs(ns::DataFrame,wfs::DataFrame)
    allCmtTracts = (wfs.ORG |> unique) ∪ (wfs.DST |> unique)
     #these tracts are not present in FluTE's `us-tracts.dat`
    notPresent = setdiff(allCmtTracts,ns.Name)
    retwfs = filter(r -> r.ORG ∉ notPresent && r.DST ∉ notPresent,eachrow(wfs))
    return retwfs
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

#=
NB! now a method for partition-to-partition commute
`ns` MUST have [:Name]; `pns` MUST have [:Name]; `prt` :Name => [ns_row_indices]
partition is reversed by a helper function, hence the need for `ns` and `pns` to have :Name
=#

function mkCmtMx(ns::DataFrame,pns::DataFrame,prt::Dict,wfs::DataFrame)
    if "Name" ∉ names(ns) || !(["ORG","DST","CMT"] ⊆ names(wfs))
        error("Wrong DF! Terminating.")
    end
    dim = nrow(pns) 
    ixs = zip(pns.Name, 1:dim) |> Dict # Name => index
   # xis = zip(1:dim, pns.Name) |> Dict # index => Name
    revp = Oboe.revexplPart(prt,ns,ixs) #generate the reverse-exploded partition
    outM = fill(0.0,(dim,dim)) # part-to-part commute matrix, init to 0.0
    for trip ∈ eachrow(wfs)
        if revp[trip.ORG] ≠ revp[trip.DST] # assuming the commute was not reflexive,
            outM[revp[trip.ORG],revp[trip.DST]]+= trip.CMT #count it in outM[pfrom,pto]
        end
    end
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
input MUST have [:Pop,:Name,:LAT,:LNG]
input SHOULD have [:IATA_Code], if absent, inserted as "dummy"
out spec: [:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG]
=#
function ns2iv_sterile(ns::DataFrame)
    out = copy(ns)
    rename!(out, :Pop => :N_i) 
    insertcols!(out,1,:id => map(select_mkid(names(out)), eachrow(out)))
    insertcols!(out,4, :S_i => out.N_i)
    insertcols!(out,5, :I_i => zeros(Int64,nrow(out)))
    insertcols!(out,6, :R_i => zeros(Int64,nrow(out)))
    if "IATA_Code" ∉ names(out)
        insertcols!(out,2,:IATA_Code => repeat(["dummy"],nrow(out)))
    end
    select(out,[:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG])
end

#quick delegate for the most common use case
ns2iv(ns::DataFrame) = ns2iv_sterile(ns) |> infect! 


#=========AUX---DBG==================#
#= This section has the functions I don't intend to move to production, which, however,
are useful enough to be pulled from disparate notebooks

Might leave them babbling out of stdout/stderr though
=#

#median of greater-than-zero elements (for the sparse stuff)
mmed(M) = filter(x -> x>0,M[1:end]) |> median
#=input: a square matrix of numbers; might add type-checking later
output: density, with a warning if nonzero found on diagonal
=#
function talkDnsy(M)
    dim = size(M,1)
    if ! (filter(x -> x > 0.0, diag(M)) |> isempty)
        println(stderr,"Warning: nonzero diagonal elements found.")
    end 
    #nonzero entries
    nzs = filter(x -> x > 0.0, M); a= nzs |> length;
    dnsy = a/ dim^2; b = dim*(dim-1) #b is max no. nonzero elements (all 'cept diag)
    println("We've got ",a," ≠0 connections")
    println("With at most ",b, " ≠0 entries, we get the density ",dnsy)
    print("min ", nzs |> minimum, "   median ", nzs |> median, "   max ", nzs |> maximum)
    return(nnonzero = a, maxnonzero= b, density = dnsy, dim = size(M,1))
end

end #end module Oboe


#========BIT=====BUCKET===========#

#collect non-equal (in view of eps) elements of 2 matrices, with their indices too
#[ (i,j,Apap[i,j] => wat.M[i,j]) for i ∈ 1:size(Apap,1) for j ∈ 1:size(Apap,1) 
#                                               if (Apap[i,j] ≉ wat.M[i,j])]
