"""
Author: Yaroslav Salii, 2020+

Main file for Oboe. Putting all the processing together,
debugging, etc.

Eventually, there'll be a couple debug scenarios and some command-line args parsing

Run me by typing `julia oboe-main-jl` at the terminal 
or `include("oboe-main.jl")` at the julia REPL

oboe-main.jl
2020-12-10  v.0.0 A Hello, World!
2021-02-02  v.0.1 Just the AP-AP to node-node processing
2021-02-19  v.0.5 Full processing, commented out; still no IO
"""

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
import Unicode

println(Oboe.callsign)

#lame logging thing
myshow = obj -> println(first(obj,5))

#initialize CLI argument parsing
s = ArgParseSettings()

@add_arg_table! s begin
    "--fips"
        nargs = '+'
        help = "FIPS codes for the states that need to be processed, or ALL for all of contiguous US"
        arg_type = AbstractString
        required = true
    "--agg"
        nargs = 1
        help = "Aggregation method {tra|cty|ap|ste}"
        arg_type = Symbol
        required = true
    "--name"
        nargs = 1
        help = "A name to be used for the output. Must not contain numbers"
        arg_type = AbstractString
        required = true
end

#Parse and process args
parsed_args = parse_args(ARGS, s)
#Get the target FIPS codes from the args, or use all states is ALL is given instead
fips = parsed_args["fips"]
if fips[1] == "ALL" && length(fips) == 1
    fips = fipsAll
end
#Get the aggregation mode and make sure it is one of the correct ones
agg = parsed_args["agg"][1]
if agg ∉ [:tra, :cty, :ap, :ste]
    throw(DomainError(agg, "The argument for agg must be {tra|cty|ap|ste}"))
end
#Finally, get the name, making sure it is upper case and has no numbers
name = Unicode.uppercase(parsed_args["name"][1])
for ch in name
    if '0' <= ch <= '9'
        throw(DomainError(name, "The name must not contain numbers"))
    end
end


iFluteFile="a~NW" * Oboe.fn.sep * Oboe.fn.fltInitSuff 
#all AP-AP travel, as list o'pairs [:ORG,:DST,:PSG], :ORG and :DST are :IATA_Code
pBTS = Oboe.grpBTS() 
#cache the smallest reasonable APs, [:IATA_Code,:LAT,:LNG,:IN.+:OUT ≥ 2500] 
APs = Oboe.censorAggFlows() # 
#read the Northwest (OR,WA) census tracts
nsRaw= Oboe.rdFluteTract(iFluteFile)
nsRaw |> myshow
#to each node, assign a designated AP from `APs` -> +[:IATA_Code]
ns = Oboe.assignDsgAPs(nsRaw,APs) 
#now find the :Pop of each APs' catchment area, and chuck that into a `Dict`
d = Oboe.mkAP_pop_dict(ns) 
#go on, compute the nodes' passenger shares (add the :shr col), with `d` in mind
ns2 = Oboe.assignPsgShares(ns,d)
ns2 |> myshow

fipsNW=["41","53"]

cmt = Oboe.rdTidyWfsByFIPS(fipsNW)
cmt |> myshow
bycty = Oboe.aggByCty(ns2)
# bycty |> myshow
byap = Oboe.aggByAP(ns2)
byste = Oboe.aggBySte(ns2)
# byste |> show
@time pNWs = Oboe.partBySte(ns2,byste)
@time pNWap = Oboe.partByAP(ns2,byap)
@time pNWc = Oboe.partByCty(ns2,bycty)


@time Act = Oboe.mkCmtMx(ns2,cmt)
# Acc = Oboe.mkCmtMx(ns2,bycty,pNWc,cmt)
# Acap= Oboe.mkCmtMx(ns2,byap,pNWap,cmt)
# Acs = Oboe.mkCmtMx(ns2,byste,pNWs,cmt)

# Apt = Oboe.mkPsgMx(ns2)
# Apc = Oboe.mkPsgMx(ns2,bycty,pNWc)
# Apap = Oboe.mkPsgMx(ns2,byap,pNWap)
# Aps = Oboe.mkPsgMx(ns2,byste, pNWs)

# inames = ["a~NW~tra","a~NW~cty","a~NW~ste"]
# inames |> println
# ivs=map(Oboe.ns2iv,[ns2,bycty,byste])
# trvs=[Apt + Act, Apc + Acc, Aps + Acs]

# for n in 1:3
#     Oboe.writeMe(inames[n],ivs[n],trvs[n])
# end
