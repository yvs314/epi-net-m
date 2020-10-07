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
=#

#= makes each `using` call re-compile this module, 
which is what I want during development; otherwise,
it continues does not reflect the changes unless the kernel is restarted =#
#__precompile__(false) 


#module, not include, to prevent multiple inculdes (oh hi #ifndef)
module Oboe 


#reading ingress $name-tract.dat, $name-wf.dat$
#reading ingress; possibly, output too
using CSV
#transforming the data in tabular form; 
using DataFrames
#read-made `mean` function
using Statistics


const callsign="This is Oboe v.0.3"
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
    dfRaw=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    names!(dfRaw, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
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

end #end module Oboe
