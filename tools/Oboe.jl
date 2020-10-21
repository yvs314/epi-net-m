#=
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
=#


#module, not include, to prevent multiple inculdes (oh hi #ifndef)
module Oboe


#reading ingress $name-tract.dat, $name-wf.dat$
#reading ingress; possibly, output too
using CSV
#transforming the data in tabular form;
using DataFrames
#read-made `mean` function
using Statistics


const callsign="This is Oboe v.0.5"
#println(callsign)

#=
1. all paths are set in view of running from epi-net-m/tools
2. data is meant to live in epi-net-m/data
=#

#all you need to know about input and output file names
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

#show the FluTE's tract filenames found in ins.ifDir, default to fn
function lsTracts(ins::NamingSpec = fn)
    filter(s::String -> endswith(s,ins.fltInitSuff),readdir(ins.ifDir))
end

#= load a FluTE census tract info into a DataFrame,
setting the types and colnames =#
function readFluteTract(ifName::String=lsTracts()[3],ins::NamingSpec=fn)
    idf=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    names!(idf, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
end

#= testing by-state aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggBySte(idf)
    by(idf,[:Ste]) do bySte
#define new rows through *named tuples*; preserves the types!
        (Pop=sum(bySte.Pop), LAT=mean(bySte.LAT), LNG=mean(bySte.LNG))
    end
end

#= testing by-county aggregation,
with dumb Euclidean centroid for geographical coordinates
Input: a FluTE $name-tracts.dat, a la [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]
=#
function aggByCty(idf)
    by(idf,[:Ste,:Cty]) do byCty
#define new rows through *named tuples*; preserves the types!
        (Pop=sum(byCty.Pop), LAT=mean(byCty.LAT), LNG=mean(byCty.LNG))
    end
end

#=======WORKING===WITH===AIRPORTS====================#

#-------TIDY---BTS---INPUT-------------------------#
#locating input files
const APdir= joinpath("..","data","by-tract","air")::String
#raw BTS data, with separate per-carrier flights
global const ifBTS=joinpath(APdir,"2019 BTS domestic.csv")::String
#raw OpenFlights AP data
global const ifAPs=joinpath(APdir,"Openflights airports.dat")::String
#TODO: consider wgetting from the original https://raw.githubusercontent.com/jpatokal/openflights/master/data/airports.dat

#=read the BTS file, and retain only :1 Passengers, :5 ORIGIN, and :7 DEST
#set the columns to [:ORG,:DST,:PSG] for uniformity=#
function rdBTS(ifName=ifBTS::String)
    rawBTS = select(CSV.read(ifName) |> DataFrame, :1,:5,:7)
    if map( x -> floor(x)==x, rawBTS[:,1]) |> all #if :PSG are Integer
        rawBTS.PASSENGERS=convert(Array{Int64,1},rawBTS.PASSENGERS)
    else #∃ non-integer passengers-per-year entry
        error("CAVEAT: *non-integer* PASSENGERS in Flights input.")
    end
    names!(rawBTS,[:PSG,:ORG,:DST])
    return rawBTS[[:ORG,:DST,:PSG]] #
end

#----------------------------------------------------#

#-------TIDY---OPENFLIGHTS---INPUT-------------------#

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

#=
pick just the `givenAPs` (by IATA_Code) from all in dfAPinfo
=#
function pickAPs(myAPs=mkFlightInfo().givenAPs::DataFrame
                , dfAPinfo=rdAPs()::DataFrame)
    out= join(myAPs, dfAPinfo,on=:IATA_Code, kind=:inner)
end

#-----BTS---AGGREGATION---ETC---------------------------#

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
#row-wise total sans disregard reflexive, `missing` if the departures are only reflexive
    ,OUT=[ (sum(M[i,:]) == M[i,i]) ? missing : (sum(M[i,:]) - M[i,i]) for i ∈ 1:size(M)[1] ]
#no. reflexive travelers, or `missing`
    ,TOUR=[Mraw[i,i] for i ∈ 1:size(Mraw)[1]] )#just the reflexive travelers
    #make the :FLOW column as :IN + :OUT; with absorbing `missing`
end

#= returns `givenAPs`, a 1-col DataFrame with all AP codes present in input
to be inner-joined ⋂ with OpenFlights to get their coordinates =#
function mkFlightInfo(idf=grpBTS()::DataFrame)
#call the auxiliary function to get the `givenAPs` and flow matrix dfA
tmp = mkMissingPairs(idf)
#add `mRts` to `idf` and transform list-of-pairs `idf` into an “adj.mx”
    dfA=unstack(vcat(idf,tmp.mRts) |> sort,:ORG,:DST,:PSG)
    #sanity check: :ORGs _[:,1] and :DSTs names(_)[2:end] are equal as sequences
    if dfA[:,1] != map(string, names(dfA))[2:end]
        error("Origin and destination names mismatch. Terminating.\n")
    end
    return(givenAPs=tmp.givenAPs,dfFlow=dfA)
end

end #end module Oboe
