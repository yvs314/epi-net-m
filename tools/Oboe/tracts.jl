"""
This module is responsible for preprocessing tract data, i.e.
things like filtering tracts and linking them with corresponding 
airport/commute data. 
"""
module Tracts
using DataFrames
#= Get Designated AP
for a node-row [:ID(:Ste,:Cty,:Tra),:Pop,:LAT,:LNG],
return the *nearest* AP's :IATA_Code and the distance :dst to it in km   

NB: worked slow-ish in Jupyter Notebook; might rework without sorting all this thing.
UPD: this slow-down was probably due to lack of caching with arguments by default
=#
function getDsgAP(node=aggBySte()[1,:]::DataFrameRow,APs=censorAggFlows()::DataFrame)
    alle=[(dst = myDist((node.LAT,node.LNG)
            ,(ap.LAT,ap.LNG))
            ,IATA_Code = ap.IATA_Code
            ,psg= ap.IN+ap.OUT)  for ap ∈ eachrow(APs)] 
    sort!(alle, by = first) # sort by distance to APs
    return (dsg=alle[1],all=alle)   
end
 
#= Assign a designated AP to each node in `nodes`, return a DF with info on both nodes and
their designated APs, [:ID,:Pop,:LAT,:LNG,:IATA_Code]; :IATA_Code for the designated AP
dbg columns: :dst distance to the dsg AP, :psg = :IN + :OUT passengers through dsg_AP
=#
function assignDsgAPs(nodes=aggBySte()::DataFrame,APs=censorAggFlows()::DataFrame)
    dsgAPs = map(n -> Oboe.getDsgAP(n,APs).dsg,eachrow(nodes)) |> DataFrame
    hcat(nodes,dsgAPs)
end

#----NODE-TO-NODE---PASSENGER--FLOWS----FROM--AP--DAILY--ENPLANEMENTS-----#
#CAVEAT: all nodes in this section must have a designated AP (:IATA_Code column)

#------------AUX-----------#
#say how many people fly through each AP; default to per-state aggregation
#`nodes` must have designated APs (:IATA_Code) and populations (:Pop)
function mkClusterPops(nodes=assignDsgAPs(aggBySte())::DataFrame)
    #group by the same designated AP: each group is this AP's __catchment area__
    gd = groupby(nodes,[:IATA_Code]) 
    #= the prime here is Σpopulation by designated AP, 
    but it's also a clusterization so here's `mean` :LAT and :LNG as dumb centroids =#
    out = combine(gd,:Pop => sum, :LAT => mean, :LNG => mean,renamecols = false)
end

#TODO: consider leaving it *inside* `assignDsgAPs`, it has to be called anyway
#make a `Dict` mapping each AP's ID (:IATA_Code) to the :Pop of its *catchment area*
function mkAP_pop_dict(nodes=assignDsgAPs(aggBySte())::DataFrame )
    cache=mkClusterPops(nodes)
    return Dict(cache.IATA_Code .=> cache.Pop)
end

#GLEaM-like aggregation: lump together all nodes in AP's catchment area (start with `census tracts`)
#difference with GLEaM: some airports may get lost, whereas GLEaM starts with APs (gotta re-check that!)
function aggByAPs()
    return assignDsgAPs(rdFluteTract()) |> mkClusterPops
end


#for a node `n`, the fraction of pop in `dsg_n`'s catchment area
#the node MUST have an :IATA_Code column (its *designated AP*)
function nodePsgShare(node = assignDsgAPs()[1,:]::DataFrameRow,dAP_pop=mkAP_pop_dict()::Dict)
    return node.Pop / dAP_pop[node.IATA_Code]
end

#apply the above to all nodes. 
#each node MUST have an :IATA_Code column (its *designated AP*)
function assignPsgShares(nodes=assignDsgAPs(aggBySte())::DataFrame,dAP_pop=mkAP_pop_dict()::Dict)
    psgShares = map(n -> Oboe.nodePsgShare(n,dAP_pop), eachrow(nodes)) #compute the shares
    out = hcat(nodes,DataFrame("shr" => psgShares)) #add them as a :shr column
end
end