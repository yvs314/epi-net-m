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
2021-07-19  v.1.0: fixed issue where unmatched tract IDs led to a crash
2021-08-17  v.1.1: updated to support the new API of Oboe.jl and FromFile module import
"""

module OboeMain

#The FIPS codes for each of the contiguous US states (and DC)
fipsAll =  ["01", "04", "05", "06", "08", "09", "10", "11", "12", "13", "16", "17", "18",
            "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
            "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "44", "45",
            "46", "47", "48", "49", "50", "51", "53", "54", "55", "56"]


using FromFile
using ArgParse
@from "./Oboe/Oboe.jl" import Oboe

export processOboe

println(Oboe.callsign)

#lame logging thing
myshow = obj -> println(first(obj,5))

function processOboe(name, agg; fips=fipsAll, useNW=false, force=false)
    #construct the output file name and make sure the files don't already exist
    iname = name * agg
    if !force && Oboe.checkIfFilesExist(iname)
        error("Output data for $iname already exists, use the force flag to overwrite")
    end

    if useNW
        iFluteFile="a~NW" * Oboe.fn.sep * Oboe.fn.fltInitSuff
        #read the census tracts corresponding to the name provided in the args
        nsRaw = Oboe.rdFluteTract(iFluteFile) 
    else
        wholeUS = Oboe.rdWholeUS()
        nsRaw = Oboe.censorFluteTractByFIPS(wholeUS, fips)
    end
    nsRaw |> myshow

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
    ns2 |> myshow

    #generate matrix of AP-to-AP flights with aux. arrays mapping indices to IATA codes
    fmx = Oboe.mkFlightMx2(pBTS, ns2.IATA_Code |> unique; daily=true)
    #generate commute table
    cmt = Oboe.rdTidyWfsByFIPS(fips, ns2)
    cmt |> myshow

    #do the processing, aggregating by tract/state/county/ap as selected in the args
    skipagg = false #this will be true if by-tract aggregation is requested
                    #and thus no additional steps are needed
    if agg == "tra"
        skipagg = true
    elseif agg == "cty"
        aggregated = Oboe.aggByCty(ns2)
        @time partitioned = Oboe.partByCty(ns2, aggregated)
    elseif agg == "ap"
        aggregated = Oboe.aggByAP(ns2)
        @time partitioned = Oboe.partByAP(ns2, aggregated)
    elseif agg == "ste"
        aggregated = Oboe.aggBySte(ns2)
        @time partitioned = Oboe.partBySte(ns2, aggregated)
    end

    @time cmtMx = skipagg ? Oboe.mkCmtMx(ns2, cmt) : Oboe.mkCmtMx(ns2, aggregated, partitioned, cmt)
    @time psgMx = skipagg ? Oboe.mkPsgMx(ns2, fmx) : Oboe.mkPsgMx(ns2, fmx, aggregated, partitioned)
    iv          = skipagg ? Oboe.ns2iv(ns2)        : Oboe.ns2iv(aggregated)


    
    iname |> println
    Oboe.writeMe(iname,iv,psgMx + cmtMx)
end

end