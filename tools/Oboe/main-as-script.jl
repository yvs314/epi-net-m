"""
Author: Yaroslav Salii, 2021+

This is a script-like version of the Oboe data processing routine.
It to be executed line-by-line from VS Code, with Julia-in-VS-Code extension

main-as-script.jl

2021-11-23 v.0.0 Collected all important things
"""

using FromFile #glue package, used to collect project source code
using DataFrames
#using CSV, Statistics 

@__DIR__

#load the project source
@from "Oboe.jl" import Oboe

# thisPath = splitpath(@__DIR__)
# projRoot = thisPath[1:findfirst(isequal("epi-net-m"), thisPath)]
# const figDir = joinpath(projRoot..., "fig")

"""
Load and pre-process the airport data
"""
rawBTS = Oboe.rdBTS()

#read Openflights.org airport data: [:IATA_Code,:LAT,:LNG,:Name]
rawAPs = Oboe.rdAPs()

#all AP-AP travel, as list o'pairs [:ORG,:DST,:PSG], :ORG and :DST are :IATA_Code
pBTS = Oboe.grpBTS(rawBTS)
#cache the smallest reasonable APs, [:IATA_Code,:LAT,:LNG,:IN.+:OUT ≥ 2500] 
APs = Oboe.getProcessedAPs(pBTS, rawAPs)

"""
Set the working area by specifying the U.S. state FIPS codes.
Only the census tracts from these states will be considered
"""
fips  = ["41","53"]
#load the whole tract database
wholeUS = Oboe.rdWholeUS()
#retain only those with matching state FIPS
nsRaw = Oboe.censorFluteTractByFIPS(wholeUS, fips)

"""
Set the aggregation mode code.
- "tra" by census tract (degenerate aggregation)
- "cty" by county
- "ap" by airport service area
- "ste" by state
"""
agg = "cty"

"""
Assign the designated airports to census tracts 
and compute the tracts' airport shares
"""
#to each node, assign a designated AP from `APs` -> +[:IATA_Code]
ns = Oboe.assignDsgAPs(nsRaw,APs)
#now find the :Pop of each APs' catchment area, and chuck that into a `Dict`
d = Oboe.mkAP_pop_dict(ns) 
#go on, compute the nodes' passenger shares (add the :shr col), with `d` in mind
ns2 = Oboe.assignPsgShares(ns,d)
#ns2 |> myshow

"""
Prepare the raw data 
- on flights (AP-to-AP flights for the designated APs)
- on commute (tract-to-tract commute for selected census tracts)
"""
#generate matrix of AP-to-AP flights with aux. arrays mapping indices to IATA codes
fmx = Oboe.mkFlightMx2(pBTS, ns2.IATA_Code |> unique; daily=true)
#generate commute table
cmt = Oboe.rdTidyWfsByFIPS(fips, ns2)
#cmt |> myshow

#do the processing, aggregating by tract/state/county/ap as selected in the args


println("Processing started. Aggregation mode [$agg] selected.")
aggnpart = Oboe.dispatchAggPart(agg) #dispatch the aggregate/partition function pair
println("Aggregating")
@time aggregated = ns2 |> aggnpart[1] #aggregate as selected by dispatcher
println("Partitioning")
@time partitioned = aggnpart[2](ns2,aggregated) #partition as selected by dispatcher

skipagg = (agg == "tra") ? true : false #no aggregation necessary if "tra" requested
println("Computing aggregated commute matrix")
@time cmtMx = skipagg ? Oboe.mkCmtMx(ns2, cmt) : Oboe.mkCmtMx(ns2, aggregated, partitioned, cmt)
println("Computing aggregated air travel matrix")
@time psgMx = skipagg ? Oboe.mkPsgMx(ns2, fmx) : Oboe.mkPsgMx(ns2, fmx, aggregated, partitioned)

#infect2! is aware of all aggregation types
iv  = Oboe.infect2!(Oboe.ns2iv_sterile(aggregated),aggType=agg) 


"""
Fix the output name
"""
# NAMEagg = "TEST" * agg

# println("Generated instance name: ",NAMEagg,"_",nrow(iv))
# println("Writing -init.csv and -trav.dat")
# Oboe.writeMe(NAMEagg,iv,psgMx + cmtMx)
# println("All done. Terminating.")