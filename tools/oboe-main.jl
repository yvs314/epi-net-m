"""
Author: Yaroslav Salii, 2020+

Main file for Oboe. Putting all the processing together,
debugging, etc.

Eventually, there'll be a couple debug scenarios and some command-line args parsing

Run me by typing `julia oboe-main-jl` at the terminal 
or `include("oboe-main.jl")` at the julia REPL

oboe-main.jl
2020-12-10  v.0.0 A Hello, World!

"""

#here's a crutch to load from current path;
#TODO: make a package
push!(LOAD_PATH,pwd())
using Oboe

println(callsign)

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
#finally, compute NODE-NODE daily air passengers
Ap =Oboe.mkPsgMx(ns2)
Oboe.talkDnsy(Ap)
