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

#here's a crutch to load from current path;
#TODO: make a package
push!(LOAD_PATH,pwd())
using Oboe

println(Oboe.callsign)

#lame logging thing
myshow = obj -> println(first(obj,5))


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
# bycty = Oboe.aggByCty(ns2)
# bycty |> myshow
# byste = Oboe.aggBySte(ns2)
# byste |> show
# pNWs = Oboe.partBySte(ns2,byste)
# pNWc = Oboe.partByCty(ns2,bycty)

# Act = Oboe.mkCmtMx(ns2,cmt)
# Acc = Oboe.mkCmtMx(ns2,bycty,pNWc,cmt)
# Acs = Oboe.mkCmtMx(ns2,byste,pNWs,cmt)

# Apt = Oboe.mkPsgMx(ns2)
# Apc = Oboe.mkPsgMx(ns2,bycty,pNWc)
# Aps = Oboe.mkPsgMx(ns2,byste, pNWs)

# inames = ["a~NW~tra","a~NW~cty","a~NW~ste"]
# inames |> println
# ivs=map(Oboe.ns2iv,[ns2,bycty,byste])
# trvs=[Apt + Act, Apc + Acc, Aps + Acs]

# for n in 1:3
#     Oboe.writeMe(inames[n],ivs[n],trvs[n])
# end
