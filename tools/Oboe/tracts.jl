"""
Submodule responsible for preprocessing tract data, i.e.
linking tracts with corresponding flight data. 
"""
module Tracts
using DataFrames
using Statistics: mean
using Distances: haversine

export assignDsgAPs,
       assignPsgShares,
       mkAP_pop_dict

#===ASSIGNING DESIGNATED AIRPORTS===#

#= Distance
Starting “haversine”/orthodromic/big-circle distance
Should get the almost true (ball vs. geoid)
when coupled with Earth's radius, or rather, AP's and loc's altitudes
CAVEAT: getting both altitudes is almost as bad 
as scrubbing “road distances” from the net
TODO: route through JuliaGeo/Geodesy
=#

"""Earth's radius in km, IAU 2015"""
const R_Earth = 6371 

"""
Calculate big-circle distance between two points on Earth.
`x` and `y` are `Union{AbstractVector{T}, NTuple{2, T}}` where T<:Real
"""
myDist(x,y) = haversine(x,y,R_Earth)

"""
Get the designated AP for a node-row `[:ID(:Ste,:Cty,:Tra),:Pop,:LAT,:LNG],
return the *nearest* AP's :IATA_Code and the distance :dst to it in km   
"""
#NB: worked slow-ish in Jupyter Notebook; might rework without sorting all this thing.
#UPD: this slow-down was probably due to lack of caching with arguments by default
function getDsgAP(node::DataFrameRow,APs::DataFrame)
    alle=[(dst = myDist((node.LAT,node.LNG)
            ,(ap.LAT,ap.LNG))
            ,IATA_Code = ap.IATA_Code
            ,psg= ap.IN+ap.OUT)  for ap ∈ eachrow(APs)] 
    sort!(alle, by = first) # sort by distance to APs
    return (dsg=alle[1],all=alle)   
end
 
"""
Assign a designated AP to each node in `nodes`, return a DF with info on both nodes and
their designated APs, `[:ID,:Pop,:LAT,:LNG,:IATA_Code]`; `:IATA_Code` for the designated AP
dbg columns: `:dst` distance to the dsg AP, `:psg` = `:IN` + `:OUT` passengers through `dsg_AP`
"""
function assignDsgAPs(nodes::DataFrame,APs::DataFrame)
    dsgAPs = map(n -> getDsgAP(n,APs).dsg,eachrow(nodes)) |> DataFrame
    hcat(nodes,dsgAPs)
end

#------------AUX-----------#
"""
Calculate how many people fly through each AP; default to per-state aggregation
`nodes` must have designated APs (:IATA_Code) and populations (:Pop)
"""
function mkClusterPops(nodes::DataFrame)
    #group by the same designated AP: each group is this AP's __catchment area__
    gd = groupby(nodes,[:IATA_Code]) 
    #= the prime here is Σpopulation by designated AP, 
    but it's also a clusterization so here's `mean` :LAT and :LNG as dumb centroids =#
    out = combine(gd,:Pop => sum, :LAT => mean, :LNG => mean,renamecols = false)
end

#TODO: consider leaving it *inside* `assignDsgAPs`, it has to be called anyway
"""Make a `Dict` mapping each AP's ID (:IATA_Code) to the :Pop of its *catchment area*"""
function mkAP_pop_dict(nodes::DataFrame)
    cache=mkClusterPops(nodes)
    return Dict(cache.IATA_Code .=> cache.Pop)
end

"""
GLEaM-like aggregation: lump together all nodes in AP's catchment area (start with *census tracts*)
difference with GLEaM: some airports may get lost, whereas GLEaM starts with APs (gotta re-check that!)
"""
function aggByAPs(tracts::DataFrame, APs::DataFrame)
    return assignDsgAPs(tracts, APs) |> mkClusterPops
end

"""
For a node `n`, find the fraction of pop in `dsg_n`'s catchment area.
#the node MUST have an `:IATA_Code` column (its *designated AP*)
"""
function nodePsgShare(node::DataFrameRow,dAP_pop::Dict)
    return node.Pop / dAP_pop[node.IATA_Code]
end

"""
Apply `nodePsgShare` to all nodes in `nodes` and store them in the `shr` column of the output. 
each node MUST have an `:IATA_Code` column (its *designated AP*)
"""
function assignPsgShares(nodes::DataFrame,dAP_pop::Dict)
    psgShares = map(n -> nodePsgShare(n,dAP_pop), eachrow(nodes)) #compute the shares
    out = hcat(nodes,DataFrame("shr" => psgShares)) #add them as a :shr column
end

end