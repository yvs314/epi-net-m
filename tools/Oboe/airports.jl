"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

Submodule responsible for preprocessing airport and flight data.

airports.jl
2021-08-13 v.0.1: First modular version
2021-08-17 v.0.2: Removed dependency on io.jl, created new func getProcessedAPs
                  to make the control flow explicit rather than rely on default arguments
2021-08-19 v.0.3: Made min_annual_boardings an optional kwarg for getProcessedAPs 
"""
module Airports
using DataFrames
using FromFile

export FlightMx,
       grpBTS,
       mkFlightMx2,
       dropSmallAPs,
       getProcessedAPs

"""The minimum number of annual boardings for airports to be counted."""
const defaultMinAnnualBoardings = 2500

"""
Data type that describes passenger flows between all airports
in a set.
"""
struct FlightMx
    """An adjacency matrix, M[i1][i2] is the flow from i1 to i2"""
    M::Matrix{Float64}
    """Maps IATA codes to indices in M"""
    ix::Dict{String, Int}
    """Maps indices in M to IATA codes"""
    xi::Dict{Int, String}
    """Lists every IATA code in the set"""
    apCodes::Vector{String}
end

#===BTS FLIGHT DATA PROCESSING===#
"""
Sum the passengers on the flights with same (ORG,DST) pairs,
sort lexicographically in (ORG,DST) order
and output the air travel graph as a *list of edges*.
"""
function grpBTS(idf::DataFrame)
    gd = groupby(idf,[:ORG,:DST])
    out  = combine(gd,:PSG => sum, renamecols = false)
    sort!(out,[:ORG,:DST]) #sort by departure/arrival AP names
end

"""
Given `apCodes [:IATA_Code]`,
outputs a `DataFrame [:IATA_Code,:IN,:OUT,:TOUR]`,
with :IN=Σ_incoming PSG, :OUT=Σ_outgoing,:TOUR=Σ_(:ORG=:DST)
in :IN and :OUT sums, `missing` is non-absorbing and the diagonal is omitted.
(opt) :TTL=:IN+:OUT.
"""
function mkAggFlows2(apCodes::Array{String}, M::Matrix{Float64})
    out=DataFrame(IATA_Code=apCodes
#col-wise total sans the reflexive, `missing` if the arrivals are only reflexive
    ,IN=[ (sum(M[:,j]) == M[j,j]) ? missing : (sum(M[:,j]) - M[j,j]) for j ∈ 1:size(M)[2] ]
#row-wise total sans the reflexive, `missing` if the departures are only reflexive
    ,OUT=[ (sum(M[i,:]) == M[i,i]) ? missing : (sum(M[i,:]) - M[i,i]) for i ∈ 1:size(M)[1] ]
#no. reflexive travelers, or `missing`
    ,TOUR=[M[i,i] for i ∈ 1:size(M)[1]] )#just the reflexive travelers
    return out
end

#===CONSTRUCTING THE FLIGHT MATRIX===#
"""
Transform a *list of arcs* into *adjacency matrix*.
input: DF must have [:ORG,:DST,:PSG] cols
"""
function mkFlightMx2(fs::DataFrame; daily=false::Bool,babble=false::Bool, init_to=0.0::Number)
    #make a sorted list of ALL the APs in `fs`
    allAPs = sort( (fs.ORG |> unique) ∪ (fs.DST |> unique))
    #delegate to the method with EXPLICIT ground set (`iretAPs`); with retaintours = true
    return mkFlightMx2(fs, allAPs; daily=daily,babble=babble,init_to=init_to,retaintours=true)
end

"""
Transform a *list of arcs* into *adjacency matrix*,
taking an *explicit* list of vertices `iretAPs`.
input: DF must have [:ORG,:DST,:PSG] cols
"""
function mkFlightMx2(fs::DataFrame, iretAPs::Array{String}; 
    daily=true::Bool # divide by 365, the input was annualized
    ,retaintours = false::Bool #flights with ORG==DST not retained by default
    ,babble=true::Bool, init_to=0.0::Number)
    
    retAPs = sort(iretAPs) #ensure AP codes are sorted, for indexing porposes
    pret = a -> a.ORG ∈ retAPs && a.DST ∈ retAPs
    if !retaintours
        pret = a -> a.ORG ≠ a.DST &&   a.ORG ∈ retAPs && a.DST ∈ retAPs
    end
    #pretain(a)  ? (retaintours) : ()
    #retain only flights (tuples :ORG,:DST,PSG) for APs ∈ retAPs; 
    retFlows = filter(pret, eachrow(fs)) |> DataFrame
    
    if babble #report if there were *any* isolated APs
        isolatedAPs = filter( a -> a ∉ retFlows.ORG && a ∉ retFlows.DST, retAPs)
        println("Found ", isolatedAPs |> length," isolated APs: ",isolatedAPs) 
    end
    
    dim = length(retAPs) 
    ix_ = zip(retAPs, 1:dim) |> Dict #get index by name 
    xi_ = zip(1:dim,retAPs) |> Dict #get name by index 
    M_ = fill(init_to,(dim,dim)) #default is 0.0; screw `missing`
    #fill the matrix with the *known* values
    for row ∈ eachrow(retFlows)
        M_[ix_[row.ORG],ix_[row.DST]] = row.PSG  
    end

    outM = daily ? map(x -> x/365,M_) : M_ #annual-to-daily, if requested
    return FlightMx(outM, ix_, xi_, retAPs)
end

#===LINKING BTS DATA WITH OPENFLIGHTS===#
#=====BTS=<-->=OPENFLIGHTS=======#
"""
Pick just the given APs (by IATA_Code) from all in dfAPinfo and
output columns as [:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country]
"""
function pickAPs(myAPs::DataFrame, dfAPinfo::DataFrame)
    out= innerjoin(myAPs, dfAPinfo,on=:IATA_Code)
    select!(out,[:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country])
end

"""
Pick just the “clean”  APs (no `missing`) from all in dfAPinfo.
CAVEAT: doesn't retain any APs that have :IN or:OUT missing or 0;
recall that 0 was identified with `missing` in mkAggFlows()=
Also remove missings from the :IN and :OUT.
"""
function pickCleanAPs(myAPs::DataFrame, dfAPinfo::DataFrame)
    pickAPs(myAPs,dfAPinfo) |> scrubAPs
end

"""Throw out the APs with missing IN or OUT enplanements"""
scrubAPs(APs::DataFrame) = filter(row -> !ismissing(row.OUT) && !ismissing(row.IN)
,eachrow(APs)) |> DataFrame

#TBD: enplanement is just OUT, so maybe not sum it with in.
"""Only retain the APs with at least `p` annual enplanements (sum :IN and :OUT)"""
function dropSmallAPs(idf::DataFrame, p=defaultMinAnnualBoardings::Number)
    if (["IN","OUT"] ⊈ names(idf)) error("wrong DF") end
    filter(row -> row.IN+row.OUT ≥ p, eachrow(idf)) |> DataFrame
end

"""
Given flight data `pBTS` and raw airport data `rawAPs`, calculate AP in- and outflows
and get rid of the APs with missing data or less than `min_annual_boardings` annual passengers.
Input:
`pBTS`: [:ORG, :DST]
`rawAPs`: [:IATA_Code,:LAT,:LNG,:Name,:City,:Country]
Output: [:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country]
"""
function getProcessedAPs(pBTS::DataFrame, rawAPs::DataFrame; min_annual_boardings=defaultMinAnnualBoardings)
    fmx = mkFlightMx2(pBTS)
    aggFlows = mkAggFlows2(fmx.apCodes, fmx.M)
    cleanAPs = pickCleanAPs(aggFlows, rawAPs)
    dropSmallAPs(cleanAPs, min_annual_boardings)
end

end