"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

Main file for Oboe. Putting all the processing together,
debugging, etc.

Eventually, there'll be a couple debug scenarios and some command-line args parsing

Run me by typing `include("oboe-main.jl").processOboe(...)` at the julia REPL

oboe-main.jl
2020-12-10  v.0.0 A Hello, World!
2021-02-02  v.0.1 Just the AP-AP to node-node processing
2021-02-19  v.0.5 Full processing, commented out; still no IO
2021-07-08  v.0.6 Implemented IO via CLI arguments
2021-07-12  v.0.7 Moved CLI parsing to a separate file
2021-04-14  v.0.8 Added support for tract datasets other than NW and checking 
                  if the output already exists
2021-07-19  v.1.0 fixed issue where unmatched tract IDs led to a crash
2021-08-17  v.1.1 updated to support the new API of Oboe.jl and FromFile module import
2021-09-08  v.1.2 isolate aggregation/partition function pair dispatch
2021-10-13  v.1.3 automatically normalize single-digit FIPS codes without a leading zero
"""

module OboeMain

"The FIPS codes for each of the contiguous US states (and DC)"
fipsAll =  ["01", "04", "05", "06", "08", "09", "10", "11", "12", "13", "16", "17", "18",
            "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
            "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "44", "45",
            "46", "47", "48", "49", "50", "51", "53", "54", "55", "56"]

"Wash. + Ore, 2072 tracts"
fipsNW = ["41","53"]

"NW's first neighbors except Cal., 2834 tracts"
fipsNW1 = ["16","32","41","53"]

"NW's 1st and 2nd neighbors except Cal., 4830 tracts, 259 ctys"
fipsNW2 = ["04","16","30","32","41","49","53","56"]

"Cal., 7038 tracts, 58 ctys"
fipsCA = ["06"]

"West Coast, Cal. + Ore. + Wash., 9910 tracts, 133 ctys"
fipsWCT= ["06","41","53"]


using FromFile
using ArgParse
@from "./Oboe/Oboe.jl" import Oboe
using DataFrames

export processOboe

println(Oboe.callsign)

#lame logging thing
myshow = obj -> println(first(obj,5))

"blueprint for getProcessedAPs() method "
#Oboe.getProcessedAPs(Oboe.rdBTS() |> Oboe.grpBTS,Oboe.rdAPs())

function processOboe(name, agg; fips=fipsAll, useNW=false, force=false)
    #construct the output file name and make sure the files don't already exist
    iname = name * agg
    if !force && Oboe.checkIfFilesExist(iname)
        error("Output data for $iname already exists, use the force flag to overwrite")
    end
    println("Preparing the tract, air travel, and commute data...")
    fips = map(Oboe.normalizeOneDigitFIPS, fips) # convert FIPS codes such as "1" to "01"
    if useNW
        iFluteFile="a~NW" * Oboe.fn.sep * Oboe.fn.fltInitSuff
        #read the census tracts corresponding to the name provided in the args
        nsRaw = Oboe.rdFluteTract(iFluteFile) 
    else
        wholeUS = Oboe.rdWholeUS()
        nsRaw = Oboe.censorFluteTractByFIPS(wholeUS, fips)
    end
    #nsRaw |> myshow

    #read flight data from BTS
    rawBTS = Oboe.rdBTS()
    #read Openflights.org airport data: [:IATA_Code,:LAT,:LNG,:Name]
    rawAPs = Oboe.rdAPs()

    #all AP-AP travel, as list o'pairs [:ORG,:DST,:PSG], :ORG and :DST are :IATA_Code
    pBTS = Oboe.grpBTS(rawBTS)
    #cache the smallest reasonable APs, [:IATA_Code,:LAT,:LNG,:IN.+:OUT ≥ 2500] 
    APs = Oboe.getProcessedAPs(pBTS, rawAPs)
    #to each node, assign a designated AP from `APs` -> +[:IATA_Code]
    ns = Oboe.assignDsgAPs(nsRaw,APs)
    #now find the :Pop of each APs' catchment area, and chuck that into a `Dict`
    d = Oboe.mkAP_pop_dict(ns) 
    #go on, compute the nodes' passenger shares (add the :shr col), with `d` in mind
    ns2 = Oboe.assignPsgShares(ns,d)
    #ns2 |> myshow

    #generate matrix of AP-to-AP flights with aux. arrays mapping indices to IATA codes
    fmx = Oboe.mkFlightMx2(pBTS, ns2.IATA_Code |> unique; daily=true)
    #generate commute table
    cmt = Oboe.rdTidyWfsByFIPS(fips, ns2)
    #cmt |> myshow

    #do the processing, aggregating by tract/state/county/ap as selected in the args
    

    println("Processing started. Aggregation mode [$agg] selected.")
    aggnpart = dispatchAggPart(agg) #dispatch the aggregate/partition function pair
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


    println("Generated instance name: ",iname,"_",nrow(iv))
    println("Writing -init.csv and -trav.dat")
    Oboe.writeMe(iname,iv,psgMx + cmtMx)
    println("All done. Terminating.")
end

"Dispatch aggregation and partition functions by aggregation type [tra | cty | ap | ste]"
function dispatchAggPart(agg::String)
    if agg == "tra"
        return (identity,Oboe.partByTra)
    elseif agg == "cty"
        return (Oboe.aggByCty,Oboe.partByCty)
    elseif agg == "ap"
        return (Oboe.aggByAP,Oboe.partByAP)
    elseif agg == "ste"
        return (Oboe.aggBySte,Oboe.partBySte)
    else
        error("Unknown aggregation requested. Terminating.\n")
    end
end

end #end module OboeMain
