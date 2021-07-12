"""
Author: Yaroslav Salii, 2020+

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
"""

module OboeMain

#The FIPS codes for each of the contiguous US states (and DC)
fipsAll =  ["01", "04", "05", "06", "08", "09", "10", "11", "12", "13", "16", "17", "18",
            "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
            "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "44", "45",
            "46", "47", "48", "49", "50", "51", "53", "54", "55", "56"]

#here's a crutch to load from current path;
#TODO: make a package
push!(LOAD_PATH,pwd())
using Oboe
using ArgParse

println(Oboe.callsign)

#lame logging thing
myshow = obj -> println(first(obj,5))

function processOboe(name, agg; fips=fipsAll, useNW=false)
    #TODO: Modify this line so that a~NW is only used when useNW is true, 
    #      and otherwise censorFluteTractByFIPS is called
    iFluteFile="a~" * name * Oboe.fn.sep * Oboe.fn.fltInitSuff 
    #all AP-AP travel, as list o'pairs [:ORG,:DST,:PSG], :ORG and :DST are :IATA_Code
    pBTS = Oboe.grpBTS() 
    #cache the smallest reasonable APs, [:IATA_Code,:LAT,:LNG,:IN.+:OUT ≥ 2500] 
    APs = Oboe.censorAggFlows()
    #read the census tracts corresponding to the name provided in the args
    nsRaw= Oboe.rdFluteTract(iFluteFile) 
    nsRaw |> myshow
    #to each node, assign a designated AP from `APs` -> +[:IATA_Code]
    ns = Oboe.assignDsgAPs(nsRaw,APs)
    #now find the :Pop of each APs' catchment area, and chuck that into a `Dict`
    d = Oboe.mkAP_pop_dict(ns) 
    #go on, compute the nodes' passenger shares (add the :shr col), with `d` in mind
    ns2 = Oboe.assignPsgShares(ns,d)
    ns2 |> myshow

    cmt = Oboe.rdTidyWfsByFIPS(fips)
    cmt |> myshow

    #do the processing, aggregating by tract/state/county/ap as selected in the args
    cmtMx = nothing
    psgMx = nothing
    iv = nothing

    if agg == "tra"
        @time cmtMx = Oboe.mkCmtMx(ns2,cmt)
        psgMx = Oboe.mkPsgMx(ns2)
        iv = Oboe.ns2iv(ns2)

    elseif agg == "cty"
        bycty = Oboe.aggByCty(ns2)
        # bycty |> myshow
        @time partCty = Oboe.partByCty(ns2,bycty)
        
        cmtMx = Oboe.mkCmtMx(ns2,bycty,partCty,cmt)
        psgMx = Oboe.mkPsgMx(ns2,bycty,partCty)
        iv = Oboe.ns2iv(bycty)

    elseif agg == "ap"
        byap = Oboe.aggByAP(ns2)
        @time partAp = Oboe.partByAP(ns2,byap)
        cmtMx = Oboe.mkCmtMx(ns2,byap,partAp,cmt)
        psgMx = Oboe.mkPsgMx(ns2,byap,partAp)
        iv = Oboe.ns2iv(byap)

    elseif agg == "ste"
        byste = Oboe.aggBySte(ns2)
        # byste |> show
        @time partSte = Oboe.partBySte(ns2,byste)
        global cmtMx = Oboe.mkCmtMx(ns2,byste,partSte,cmt)
        global psgMx = Oboe.mkPsgMx(ns2,byste, partSte)
        global iv = Oboe.ns2iv(byste)
    end

    iname = name * agg
    iname |> println

    Oboe.writeMe(iname,iv,psgMx + cmtMx)
end

end