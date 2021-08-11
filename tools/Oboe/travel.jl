"""
This module is responsible for producing the commute and air travel matrices.
"""
module Travel
using DataFrames
using FromFile

@from "airports.jl" using Airports: mkFlightMx2, grpBTS
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

#NB! `ns`:[:ID,:Pop,:IATA_Code,:shr]
#NB! if provided, `pns` MUST have [:Name]; `prt` :Name => [ns_row_indices]
#    mkPsgMx(ns[, pns, prt][, force_recompute=true|false])
#Calculates the air passenger flow in both directions for all node pairs in `ns`.
# - If `pns` and `prt` are not provided, returns a matrix `M` such that 
#   `M[i1][i2]` is the total air travel flow from `ns[i1]` to `ns[i2]`
# - If `pns` and `prt` are provided, returns a matrix `M` such that 
#   `M[i1][i2]` is the total air travel flow from `pns[i1]` to `pns[i2]`,
#   obtained by summing the flows from and to each node in the corresponding
#   partitions of `prt`.
function mkPsgMx(ns::DataFrame,
                 pns::Union{DataFrame,Nothing}=nothing,
                 prt::Union{Dict, Nothing}=nothing;
                 force_recompute::Bool=false)
    #first and foremost get the list of all APs present in ns... 
    retAPs = ns.IATA_Code |> unique
    #...then get the daily AP-to-AP flows for the designated APs, with IATA_Code as index
    fmx = mkFlightMx2(grpBTS(),retAPs; daily=true)
    #determine if we're dealing with aggregated data
    if pns !== nothing && prt !== nothing
        aggregated = true
        #assuming force_recompute is disabled, try to detect by-AP aggregation..
        if !force_recompute && "IATA_Code" ∈ names(pns) && pns.IATA_Code == sort(retAPs)
            #...and if it is, just return the by-AP travel matrix
            println("By-AP aggregation engaged. Using AP-to-AP travel matrix directly.")
            return fmx.M
        end
        #the number of partitions is the size of our output matrix
        dim = length(prt)
        #build a dict mapping the name of a partition to its index in pns
        partIndexByName = Dict{String, Int}(
            partName => findfirst(isequal(partName), pns.Name)[1] for partName in keys(prt)
        )
        #build an array mapping each node in ns to the index of its partition in pns
        partIndexByNode = zeros(Int, nrow(ns))
        #prt contains the ns indices of each node in each partition
        #so we loop through prt assigning the corresponding part. index to each node index
        for partName ∈ keys(prt), nodeIndex ∈ prt[partName]
            partIndexByNode[nodeIndex] = partIndexByName[partName]
        end
    else
        aggregated = false
        #if we're dealing with non-aggregated tract data,
        #this is treated as a trivial case where each node is its own partition
        dim = nrow(ns)
        partIndexByNode = 1:dim
    end
    #attach partition indices to ns so that we remember where to save results
    insertcols!(ns, :partindex => partIndexByNode)
    #group the nodes into subDFs by their designated airport
    groupedbyAP = groupby(ns, :IATA_Code)
    # numAPs = length(groupedbyAP)
    #from this, get list of IATA codes of unique airports that have nodes assigned
    # retAPs = ns.IATA_Code |> unique
    #initialize output matrix
    outM = fill(0.0, (dim,dim))
    computePsgFlows!(outM, retAPs, fmx, groupedbyAP, accumulate=aggregated)
    select!(ns, Not(:partindex)) #drop the partindex column to keep ns clean
    return outM
end


#For each combination (i.e. unordered pair with distinct entries)
#of groups in `groupedByAP`, and for each node pair within these combinations,
#calculates the passenger flows for that pair and accumulates them in `outM`
#at a location determined by the `partindex` column.
#So M[i1, i2] is the total passenger flow from all nodes that have the `partindex` i1
#to all nodes that have the `partindex` i2.
#
#If `accumulate` is `false`, values in `outM` will be overwritten rather than accumulated,
#which is useful for saving on memory accesses in the trivial case where each partition
#has only a single node.
@inline function computePsgFlows!(outM::Matrix{Float64},
                                  retAPs::AbstractArray{String},
                                  fmx::NamedTuple,
                                  groupedByAP::GroupedDataFrame;
                                  accumulate::Bool=true)
    numAPs = length(retAPs)
    #for each combination of airports (AP1, AP2)...
    #(this weird iteration is to avoid repeating pairs)
    for i1 ∈ 1:numAPs, i2 ∈ (i1 + 1):numAPs
        # get IATA codes of AP1 and AP2
        iata1 = retAPs[i1]
        iata2 = retAPs[i2]
        # get AP1->AP2 and AP2->AP1 passenger numbers
        flow1to2 = fmx.M[fmx.ix[iata1], fmx.ix[iata2]]
        flow2to1 = fmx.M[fmx.ix[iata2], fmx.ix[iata1]]
        # skip over airport pairs that have no travel between them
        if flow1to2 == 0.0 && flow2to1 == 0.0
            continue
        end
        grp1 = groupedByAP[i1]
        grp2 = groupedByAP[i2]
        #pre-collect relevant columns to optimize performance
        indices1 = collect(grp1.partindex) :: Vector{Int}
        indices2 = collect(grp2.partindex) :: Vector{Int}
        shr1     = collect(grp1.shr)       :: Vector{Float64}
        shr2     = collect(grp2.shr)       :: Vector{Float64}
        #for each pair of nodes (n1 ∈ AP1, n2 ∈ AP2)...
        for n1 ∈ 1:nrow(grp1), n2 ∈ 1:nrow(grp2)
            #calculate psg flow n1->n2 and n2->n1 and store the result
            outi1 = indices1[n1]
            outi2 = indices2[n2]
            #do not count travel within the same subdivision
            if outi1 == outi2
                continue
            end
            result1to2 = psg(shr1[n1], shr2[n2], flow1to2)
            result2to1 = psg(shr2[n2], shr1[n1], flow2to1)
            #if we're accumulating, add the results to corresponding outM location...
            if accumulate
                outM[outi1, outi2] += result1to2
                outM[outi2, outi1] += result2to1
            #..otherwise overwrite
            else
                outM[outi1, outi2] = result1to2
                outM[outi2, outi1] = result2to1
            end
        end
    end
end

function psg(fromshr, toshr, totalflow)
    fromshr * toshr * totalflow
end


end