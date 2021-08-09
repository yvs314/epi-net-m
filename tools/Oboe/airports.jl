"""
This module is responsible for preprocessing airport and flight data.
"""
module Airports
using DataFrames

#=======WORKING===WITH===AIRPORTS====================#


#-------TIDY---OPENFLIGHTS---INPUT-------------------#

#move these consts to BASE?
const apAllColNames = [:ID,:Name,:City,:Country,:IATA_Code,:ICAO_Code,:LAT,:LNG,:Altitude,:Timezone,:Daylight_Savings,:TZ,:Type,:Source]
#reordered in the order of necessity; 3-letter IATA code as ID
const apRetainedColNames=[:IATA_Code,:LAT,:LNG,:Name,:City,:Country]

#=
read the OpenFlights.org's airports.dat, and retain only the columns
I feel I may use; put :IATA_Code,:LAT,:LNG first,
these are definitely useful
CAVEAT: some missing values in :City, some weird values in :Country
TODO: resolve this caveat
=#
function rdAPs(ifName=ifAPs::String)
    out=CSV.File(ifName,header=false) |> DataFrame
    rename!(out, apAllColNames)
    #retain only the useful columns
    select!(out,apRetainedColNames)
end


#============BTS===FLIGHT===DATA===PROCESSING============#

#---READ---BTS--FLIGHTS-------------------------#
#=read the BTS file, and retain only :1 Passengers, :5 ORIGIN, and :7 DEST
#set the columns to [:ORG,:DST,:PSG] for uniformity=#
function rdBTS(ifName=ifBTS::String)
    rawBTS = select(CSV.read(ifName,DataFrame), :1,:5,:7)
    if map( x -> floor(x)==x, rawBTS[:,1]) |> all #if :PSG are Integer
        rawBTS.PASSENGERS=convert(Array{Int64,1},rawBTS.PASSENGERS)
    else #∃ non-integer passengers-per-year entry
        error("CAVEAT: *non-integer* PASSENGERS in Flights input.")
    end
    rename!(rawBTS,[:PSG,:ORG,:DST])
    return select(rawBTS,[:ORG,:DST,:PSG]) #
end


#---BTS---TIDY--AGGREGATE---ETC---------------------------#

#=
sum the passengers on the flights with same (ORG,DST) pairs
sort lexicographically in (ORG,DST) order
final output is the air travel graph as a *list of edges*
=#
function grpBTS(idf=rdBTS()::DataFrame)
    gd = groupby(idf,[:ORG,:DST])
    out  = combine(gd,:PSG => sum, renamecols = false)
    sort!(out,[:ORG,:DST]) #sort by departure/arrival AP names
end


#=
Transform a *list of arcs* into *adjacency matrix*
input: DF must have [:ORG,:DST,:PSG] cols
=#
function mkFlightMx2(fs=grpBTS()::DataFrame; daily=false::Bool,babble=false::Bool, init_to=0.0::Number)
    #make a sorted list of ALL the APs in `fs`
    allAPs = sort( (fs.ORG |> unique) ∪ (fs.DST |> unique))
    #delegate to the method with EXPLICIT ground set (`iretAPs`); with retaintours = true
    return mkFlightMx2(fs, allAPs; daily=daily,babble=babble,init_to=init_to,retaintours=true)
end

#takes *explicit* list of vertices `iretAPs`
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
    return (M=outM, ix = ix_, xi = xi_,apCodes = retAPs)
end

#=
in: apCodes [:IATA_Code], 
output: a `df` [:IATA_Code,:IN,:OUT,:TOUR],
with :IN=Σ_incoming PSG, :OUT=Σ_outgoing,:TOUR=Σ_(:ORG=:DST)
in :IN and :OUT sums, `missing` is non-absorbing and the diagonal is omitted
(opt) :TTL=:IN+:OUT,
=#
function mkAggFlows2(apCodes::Array{String}=mkFlightMx2().apCodes, M=mkFlightMx2().M::Matrix)
        out=DataFrame(IATA_Code=apCodes
    #col-wise total sans the reflexive, `missing` if the arrivals are only reflexive
        ,IN=[ (sum(M[:,j]) == M[j,j]) ? missing : (sum(M[:,j]) - M[j,j]) for j ∈ 1:size(M)[2] ]
    #row-wise total sans the reflexive, `missing` if the departures are only reflexive
        ,OUT=[ (sum(M[i,:]) == M[i,i]) ? missing : (sum(M[i,:]) - M[i,i]) for i ∈ 1:size(M)[1] ]
    #no. reflexive travelers, or `missing`
        ,TOUR=[M[i,i] for i ∈ 1:size(M)[1]] )#just the reflexive travelers
    return out
end


#------BTS----CENSORING--------------------#

#TBD: enplanement is just OUT, so maybe not sum it with in.
#TBD: rename to dropSmallAPs
#only retain the APs with at least `p` annual enplanements (sum :IN and :OUT)
function censorAggFlows(p=2500::Number,idf=pickCleanAPs()::DataFrame)
    if (["IN","OUT"] ⊈ names(idf)) error("wrong DF") end
    filter(row -> row.IN+row.OUT ≥ p, eachrow(idf)) |> DataFrame
end

#throw out the APs with missing IN or OUT enplanements
scrubAPs(APs=mkAggFlows2()::DataFrame) = filter(row -> !ismissing(row.OUT) && !ismissing(row.IN)
,eachrow(APs)) |> DataFrame

#===========INTERCONNECT===================#

#=====BTS=<-->=OPENFLIGHTS=======#


#pick just the given APs (by IATA_Code) from all in dfAPinfo
#output columns as [:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country]
function pickAPs(myAPs=mkAggFlows2()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
    out= innerjoin(myAPs, dfAPinfo,on=:IATA_Code)
    select!(out,[:IATA_Code,:LAT,:LNG,:IN,:OUT,:TOUR,:Name,:City,:Country])
end

#= pick just the “clean”  APs (no `missing`) from all in dfAPinfo
CAVEAT: doesn't retain any APs that have :IN or:OUT missing or 0;
recall that 0 was identified with `missing` in mkAggFlows()=
Also remove missings from the :IN and :OUT
=#
function pickCleanAPs(myAPs=mkAggFlows2()::DataFrame, dfAPinfo=rdAPs()::DataFrame)
    pickAPs(myAPs,dfAPinfo) |> scrubAPs
end


end