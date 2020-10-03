#= 
Author: Yaroslav Salii, 2020.

This is Oboe-Mangle, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS

Reading the FluTE data and transforming it to my input format.
Might separate my input format definitions later.

oboe-main.jl v.0.1: "From scripts to proper code" edition
let's migrate from Jupyter NB, I say.
=#

#module, not include, to prevent multiple inculdes (oh hi #ifndef)
module Oboe 

#reading ingress $name-tract.dat, $name-wf.dat$
#reading ingress; possibly, output too
using CSV
#transforming the data in tabular form; 
using DataFrames


global const callsign="This is Oboe v.0.1"
println(callsign)

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
global const fn=NamingSpec("-","_"
    ,joinpath("..","data","by-tract","flute")
    ,joinpath("..","data","by-tract")
    ,"tracts.dat","init.csv")

end
