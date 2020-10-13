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


const callsign="This is Oboe v.0.4"
#println(callsign)

#=
1. all paths are set in view of running from epi-net-m/tools
2. data is meant to live be in epi-net-m/data
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
    # # (rel) path to BTS Flights
    # ifBTS::String
    # # (rel) path to OpenFlights airports.dat
    # ifOAP::String
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
    #testing the exit on non-integer passengers
    #rawBTS.PASSENGERS[1]=1.1
    if map( x -> floor(x)==x, rawBTS[:,1]) |> all
        rawBTS.PASSENGERS=convert(Array{Int64,1},rawBTS.PASSENGERS)
    else
        println("CAVEAT: *non-integer* PASSENGERS in Flights input.")
    end
    names!(rawBTS,[:PSG,:ORG,:DST])
    return rawBTS[[:ORG,:DST,:PSG]] #
end

#----------------------------------------------------#

#-------TIDY---OPENFLIGHTS---INPUT-------------------#

const apColNames = [:ID,:Name,:City,:Country,:IATA_Code,:ICAO_Code,:LAT,:LNG,:Altitude,:Timezone,:Daylight_Savings,:TZ,:Type,:Source]

#just read the OpenFlights.org's airports.dat, giving proper col names
function rdAPs(ifName=ifAPs::String)
    out=CSV.File(ifName,header=false) |> DataFrame
  #  ,types=[Int,String,String,String,String,String,Float64,Float64,Int,Float64,String,String,String]) |> DataFrame
    names!(out, apColNames)
    #gotta drop the unnecessary columns now
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
    sort!(out,[:ORG,:DST])
end

#= returns `givenAPs`, a 1-col DataFrame with all AP codes present in input
to be inner-joined ⋂ with OpenFlights to get their coordinates =# 
function mkFlightInfo(idf=grpBTS()::DataFrame)
    #separate Org APs and Dest APs
    uOrgs = select(idf, :ORG) |> sort |> unique 
    uDsts = select(idf, :DST) |> sort |> unique
    #give "em the same column name, to have convenient "join" ops
    rename!(uOrgs,Dict(:ORG=> :apName))
    rename!(uDsts,Dict(:DST=> :apName))
    #let's find 2-way APs and all APs
    twoWayAPs = join(uOrgs,uDsts, on = :apName, kind = :inner) #inner join: orgs ⋂ dsts
    allAPs = join(uOrgs,uDsts, on = :apName, kind = :outer) #outer join: orgs ⋃ dsts
    
    # missing origins: allAPs ∖ uOrgs; :anti-join for ∖setminus
    mOrgs = join(allAPs,uOrgs, on = :apName, kind = :anti)
    # missing dests: allAPs ∖ uDsts; :anti-join for ∖setminus
    mDsts = join(allAPs,uDsts, on = :apName, kind = :anti)
    return(givenAPs=allAPs,orgs=uOrgs,dests=uDsts)    
end

end #end module Oboe
