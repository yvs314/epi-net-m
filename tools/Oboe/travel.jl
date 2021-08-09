"""
This module is responsible for producing the commute and air travel matrices.
"""
module Travel
using DataFrames
using FromFile

@from "airports.jl" using Airports: mkFlightMx2, grpBTS
@from "tracts.jl" using Tracts: assignPsgShares
@from "aggregation.jl" using Aggregation: revexplPart

export mkCmtMx,
       mkPsgMx

#=return wfs excluding the commutes for tracts not present in `ns`
ns MUST have [:Name]; wfs MUST have [:ORG,:DST,:CMT]
=#
function scrubWfs(ns::DataFrame,wfs::DataFrame)
    allCmtTracts = (wfs.ORG |> unique) ∪ (wfs.DST |> unique)
     #these tracts are not present in FluTE's `us-tracts.dat`
    notPresent = setdiff(allCmtTracts,ns.Name)
    retwfs = filter(r -> r.ORG ∉ notPresent && r.DST ∉ notPresent,eachrow(wfs))
    return retwfs
end



#--------COMMUTER---FLOW---MATRIX--------------#


function mkCmtMx(ns::DataFrame,wfs::DataFrame)
    if "Name" ∉ names(ns) || !(["ORG","DST","CMT"] ⊆ names(wfs))
        error("Wrong DF! Terminating.")
    end
    dim = nrow(ns) 
    outM = fill(0.0,(dim,dim)) #commute matrix, node-to-node, ordered as `ns`
    ixs = zip(ns.Name, 1:dim) |> Dict #get index by name
    xis = zip(1:dim, ns.Name) |> Dict #get name by index
    for trip ∈ eachrow(wfs) #doggedly fill each trip into the matrix
        outM[ixs[trip.ORG],ixs[trip.DST]]=trip.CMT
    end #reflexive commute is removed beforehand

    return outM
end

#=
NB! now a method for partition-to-partition commute
`ns` MUST have [:Name]; `pns` MUST have [:Name]; `prt` :Name => [ns_row_indices]
partition is reversed by a helper function, hence the need for `ns` and `pns` to have :Name
=#

function mkCmtMx(ns::DataFrame,pns::DataFrame,prt::Dict,wfs::DataFrame)
    if "Name" ∉ names(ns) || !(["ORG","DST","CMT"] ⊆ names(wfs))
        error("Wrong DF! Terminating.")
    end
    dim = nrow(pns) 
    ixs = zip(pns.Name, 1:dim) |> Dict # Name => index
   # xis = zip(1:dim, pns.Name) |> Dict # index => Name
    revp = revexplPart(prt,ns,ixs) #generate the reverse-exploded partition
    outM = fill(0.0,(dim,dim)) # part-to-part commute matrix, init to 0.0
    for trip ∈ eachrow(wfs)
        if revp[trip.ORG] ≠ revp[trip.DST] # assuming the commute was not reflexive,
            outM[revp[trip.ORG],revp[trip.DST]]+= trip.CMT #count it in outM[pfrom,pto]
        end
    end
    return outM
end

#------AIR---PASSENGER---FLOW---MATRIX----------------------#
#NB! `nodes`:[:ID,:Pop,:IATA_Code,:shr]
function mkPsgMx(ns=assignPsgShares()::DataFrame)
    retAPs = ns.IATA_Code |> unique #the APs that are designated for at least one `node`
    #get the daily AP-to-AP flows for the designated APs, with IATA_Code as index
    aps = mkFlightMx2(grpBTS(),retAPs; daily=true)  
    dim = nrow(ns) #final output matrix is [dim × dim]
    outM = fill(0.0, (dim,dim))
    #for each [from,to] pair, delegate to aux function psg
    for from ∈ 1:dim, to ∈ 1:dim 
        outM[from,to] = psg(from,to,ns,aps) #NB! reflexive flights are set to 0.0
    end 
    return outM
end

#=
NB! now a method for partition-to-partition flights
`ns` MUST have [:IATA_Code,:shr]; `pns` MUST have [:Name]; `prt` :Name => [ns_row_indices]
=#
function mkPsgMx(ns::DataFrame,pns::DataFrame,prt::Dict;force_recompute=false)
    retAPs = ns.IATA_Code |> unique |> sort #the APs that are designated for at least one `node`
    aps = mkFlightMx2(grpBTS(),retAPs; daily=true)  #get the daily AP-to-AP flows for the designated APs
    dim = nrow(pns) #final output matrix is [dim × dim], for nodes in `pns`
    outM = fill(0.0, (dim,dim))
    if("IATA_Code" ∈ (pns |> names) && #detect by-AP aggregation
        !force_recompute && #explicit request not to recompute the matrix
        retAPs == pns.IATA_Code) #ensure by-AP aggregation was correct for `ns`
        println("By-AP aggregation engaged. Using AP-to-AP travel matrix directly.")
        outM = aps.M #the diagonal is zero by definition of mkFlightMx2()
    else
    #proceed column-wise
        for pto ∈ 1:dim, pfrom ∈ 1:dim
            if pfrom ≠ pto #no reflexive air travel (reflexive defaults to 0.0 by initialization)
                outM[pfrom,pto] = sum(psg(efrom,eto,ns,aps) #sum the travel between consituents
                    for eto ∈ prt[pns.Name[pto]], 
                        efrom ∈ prt[pns.Name[pfrom]] )
            end
        end
    end
    return outM
end


#---AUX::---AIR---PASSENGER---FLOW---MATRIX----------------------#
#APs MUST be a NamedTuple (M,ix,xi), as returned by mkFlightMx2
#ns MUST have [:shr,:IATA_Code]
function psg(from::Integer,to::Integer,ns::DataFrame,APs)
    if ns.IATA_Code[from] ≠ ns.IATA_Code[to]
        return ns.shr[from] * ns.shr[to] * 
            APs.M[APs.ix[ns.IATA_Code[from]],APs.ix[ns.IATA_Code[to]]] 
    else return 0.0
    end
end

end