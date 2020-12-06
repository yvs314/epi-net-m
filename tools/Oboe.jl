"""
Author: Yaroslav Salii, 2020.

This is Oboe-Mangle, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS

Reading the FluTE data and transforming it to my input format.
Might separate my input format definitions later.

oboe-main.jl v.0.1: "From scripts to proper code" edition
             v.0.2: switch to `module`
             v.0.3: got as far as reading FluTE tracts
             v.0.4: added basic by-county and by-state IV aggregation
             v.0.5: join (intersect) BTS with OpenFlights
             v.0.6: added designated AP assignment
"""


#module, not `include`, to prevent multiple inculdes (oh hi #ifndef)
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

const callsign="This is Oboe v.0.6"
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
# sample usage: Oboe.lsTracts()[4] |> Oboe.readFluteTract |> Oboe.aggBySte
# a *node* (tract) must have the fields :Pop,:LAT,:LNG

#show the FluTE's tract filenames found in ins.ifDir, default to fn
function lsTracts(ins::NamingSpec = fn)
    filter(s::String -> endswith(s,ins.fltInitSuff),readdir(ins.ifDir))
end

#= load a FluTE census tract info into a DataFrame,
setting the types and colnames 
Default to all-US tracts
=#
function readFluteTract(ifName::String=lsTracts()[4],ins::NamingSpec=fn)
    idf=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    names!(idf, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
    return idf
end

#--------AGGREGATE---TRACT-LIKE---DATA--------------#
#= testing by-state aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggBySte(idf=readFluteTract()::DataFrame)
    by(idf,[:Ste]) do bySte
#define new rows through *named tuples*; preserves the types!
        (Pop=sum(bySte.Pop), LAT=mean(bySte.LAT), LNG=mean(bySte.LNG))
    end
end

#= testing by-county aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggByCty(idf=readFluteTract()::DataFrame)
    by(idf,[:Ste,:Cty]) do byCty
#define new rows through *named tuples*; preserves the types!
        (Pop=sum(byCty.Pop), LAT=mean(byCty.LAT), LNG=mean(byCty.LNG))
    end
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
    names!(out, apAllColNames)
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
    names!(rawBTS,[:PSG,:ORG,:DST])
    return rawBTS[[:ORG,:DST,:PSG]] #
end


#---BTS---TIDY--AGGREGATE---ETC---------------------------#

#=
sum the passengers on the flights with same (ORG,DST) pairs
(moves Passengers to the 3rd column)
sort lexicographically in (ORG,DST) order
set the column names to [:ORG,:DST,:PSG]
final output is the air travel graph as a *list of edges*
=#
function grpBTS(idf=rdBTS()::DataFrame)
    out = by(idf, [:ORG,:DST]) do flights
        [sum(flights.PSG)]
    end
    names!(out, [:ORG,:DST,:PSG]) #restore :PSG's name from :x1
    sort!(out,[:ORG,:DST]) #sort by departure/arrival AP names
end

#=
Take idf::[:ORG,:DST,:PSG], return the “missing pairs” such that 
∀i: out[i,_,_] (i ∉ idf.ORG) ∧ (i ∈ idf.DST)
∀j: out[_,j,_] (j ∉ idf.DST) ∧ (i ∈ idf.ORG)
and also “all-pairs”, idf.ORG ⋃ idf.DST 
=#
function mkMissingPairs(idf::DataFrame)
    #pick :ORG and :DST APs, sort them, and set col name to :IATA_Code
    uOrgs = select(idf, :ORG) |> unique |> sort |> (df -> names!(df,[:IATA_Code]))
    uDsts = select(idf, :DST) |> unique |> sort |> (df -> names!(df,[:IATA_Code]))
#list all APs *mentioned*, whether normal, reflexive, or in/out-only
#outer join: orgs ⋃ dsts    
    allAPs = join(uOrgs,uDsts, on = :IATA_Code, kind = :outer) |> sort
# missing origins: allAPs ∖ uOrgs; :anti-join for ∖setminus
    mOrgs = join(allAPs,uOrgs, on = :IATA_Code, kind = :anti)
# missing dests: allAPs ∖ uDsts; :anti-join for ∖setminus
    mDsts = join(allAPs,uDsts, on = :IATA_Code, kind = :anti)
#make up the missing (:ORG,:DST) pairs, i.e., mOrgs × mDests
    odf = join(mOrgs,mDsts, kind = :cross, makeunique=true)
    names!(odf,[:ORG,:DST]) #restore the [:ORG,:DST] names
    #add the dummy :PSG column (all `missing`), after :ORG and :DST
    insertcols!(odf,3,:PSG => repeat([missing::Union{Int64,Missing}], nrow(odf)))
    return (mRts=odf,givenAPs=allAPs)
end



#= returns `givenAPs`, a 1-col DataFrame with all AP codes present in input
to be inner-joined ⋂ with OpenFlights to get their coordinates =#
function mkFlightInfo(idf=grpBTS()::DataFrame)
#call the auxiliary function to get the `givenAPs` and flow matrix dfA
    tmp = mkMissingPairs(idf)
#add `mRts` to `idf` and transform list-of-pairs `idf` into an “adj.mx”
    dfA=unstack(vcat(idf,tmp.mRts) |> sort, :ORG,:DST,:PSG)
    #sanity check: :ORGs _[:,1] and :DSTs names(_)[2:end] are equal as sequences
    if dfA[:,1] != map(string, names(dfA))[2:end]
        error("Origin and destination names mismatch. Terminating.\n")
    end
    return (givenAPs=tmp.givenAPs,dfFlow=dfA)
end

#= take dFlow::[:ORG,:AP1,:AP2,...,AP_n], [n × (n+1)],
a matrix-like `df`, with AP names in its 1st col, _[:,1]
return a `df` out::[:IATA_Code,:FLOW,:IN,:OUT,:TOUR],
with :FLW=:IN+:OUT,:IN=Σ_incoming PSG, :OUT=Σ_outgoing,:TOUR=Σ_(:ORG=:DST)
in :IN and :OUT sums, `missing` is non-absorbing and the diagonal is omitted
=#
function mkAggFlows(dfFlow=mkFlightInfo().dfFlow::DataFrame)
    Mraw = convert(Matrix,dfFlow[:,2:end])
    # identify `missing` with 0; set 0 everywhere
    M = map( x -> ismissing(x) ? 0 : x, Mraw)
    out=DataFrame(IATA_Code=dfFlow[:,1]
    #col-wise total sans the reflexive, `missing` if the arrivals are only reflexive
        ,IN=[ (sum(M[:,j]) == M[j,j]) ? missing : (sum(M[:,j]) - M[j,j]) for j ∈ 1:size(M)[2] ]
    #row-wise total sans the reflexive, `missing` if the departures are only reflexive
        ,OUT=[ (sum(M[i,:]) == M[i,i]) ? missing : (sum(M[i,:]) - M[i,i]) for i ∈ 1:size(M)[1] ]
    #no. reflexive travelers, or `missing`
        ,TOUR=[Mraw[i,i] for i ∈ 1:size(Mraw)[1]] )#just the reflexive travelers
        #make the :FLOW column as :IN + :OUT; with absorbing `missing`
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
scrubAPs(APs=mkAggFlows()::DataFrame) = filter(row -> !ismissing(row.OUT) && !ismissing(row.IN)
,eachrow(APs)) |> DataFrame

#===========INTERCONNECT===================#

#=====BTS=<-->=OPENFLIGHTS=======#


#pick just the given APs (by IATA_Code) from all in dfAPinfo
#output columns as [:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country]
function pickAPs(myAPs=mkAggFlows()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
    out= innerjoin(myAPs, dfAPinfo,on=:IATA_Code)
    select!(out,[:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country])
end

#= pick just the “clean”  APs (no `missing`) from all in dfAPinfo
CAVEAT: doesn't retain any APs that have :IN or:OUT missing or 0;
recall that 0 was identified with `missing` in mkAggFlows()=
Also remove missings from the :IN and :OUT
=#
function pickCleanAPs(myAPs=mkAggFlows()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
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
function getDsgAP(node=aggBySte()[1,:]::DataFrameRow,APs=pickCleanAPs()::DataFrame)
    alle=[(dst = myDist((node.LAT,node.LNG)
            ,(ap.LAT,ap.LNG))
            ,IATA_Code = ap.IATA_Code
            ,psg= ap.IN+ap.OUT)  for ap ∈ eachrow(APs)] 
    sort!(alle, by = first) # sort by distance to APs
    return (dsg=alle[1],all=alle)   
end

# Apply getDsgAP to each node, and write that into its node
function assignDsgAPs(nodes=aggBySte()::DataFrame,APs=censorAggFlows()::DataFrame)
    dsgAPs = map(n -> Oboe.getDsgAP(n,APs).dsg,eachrow(nodes)) |> DataFrame
    hcat(nodes,dsgAPs)
end

#----NODE-TO-NODE---PASSENGER--FLOWS----FROM--AP--DAILY--ENPLANEMENTS-----#
#say how many people fly through each AP; default to per-county aggregation
#`nodes` must have designated APs (:IATA_Code) and populations (:Pop)
function mkClusterPops(nodes=assignDsgAPs(aggByCty())::DataFrame)
    by(nodes,[:IATA_Code]) do byDsgAP
        #define new rows through *named tuples*; preserves the types!
                (Pop=sum(byDsgAP.Pop), LAT=mean(byDsgAP.LAT), LNG=mean(byDsgAP.LNG))
            end
      #return nodes
end

#GLEaM-like aggregation (lump together all nodes in AP's catchment area)
#difference: some airports may get lost, whereas GLEaM starts with APs (gotta re-check that!)
function aggByAPs(nodes=assignDsgAPs(readFluteTract())::DataFrame)
    by(nodes,[:IATA_Code]) do byDsgAP
        #define new rows through *named tuples*; preserves the types!
                (Pop=sum(byDsgAP.Pop), LAT=mean(byDsgAP.LAT), LNG=mean(byDsgAP.LNG))
            end
      #return nodes
end

#for a node `n`, the fraction of passengers in `dsg_n` to/from `n`
function nodePsgShare()
    return 0.0
end

end #end module Oboe
