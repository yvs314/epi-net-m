"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

Submodule responsible for grouping and aggregating tracts to various levels.

aggregation.jl
2021-08-13 v.0.1: First modular version
"""
module Aggregation
using DataFrames
using Statistics: mean

export aggByCty,
       aggBySte,
       aggByAP,
       partByCty,
       partByAP,
       partBySte

 #NB! :Name is a String
 """Return a name-generating function for partitions based on the column names `nms`"""
 function select_mkName(nms::Array{String})
    if ["Ste","Cty","Tra"] ⊆ nms #node's a census tract
         return r-> join([r.Ste,r.Cty,r.Tra],"~")
    elseif ["Ste","Cty"] ⊆ nms #node's a county
         return r-> join([r.Ste,r.Cty],"~")
    elseif  ["Ste"] ⊆ nms #node's a state
         return r-> string(r.Ste) #ensure it's output as a string
     elseif ["IATA_Code"] ⊆ nms #node is Voronoi cell around an AP
         return r -> r.IATA_Code
     else
        error("Nodes header not recognized. Can't generate a :Name without FIPS or IATA_Code.")
     end
 end

"""
By-state PRE-aggregation with dumb Euclidean centroid for geographical coordinates.
Input: a FluTE \$name-tracts.dat, a la `[:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]``
"""
#TODO: put the by--end output into a variable and add post-processing (adding the ID)
function aggBySte(idf::DataFrame;make_names=true)
    gd = groupby(idf,[:Ste]) #group by U.S. State FIPS
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
end

""" 
By-county PRE-aggregation with dumb Euclidean centroid for geographical coordinates.
Input: a FluTE \$name-tracts.dat, a la `[:Ste,:Cty,:Tra,:Pop,:LAT,:LNG]`
"""
function aggByCty(idf::DataFrame;make_names=true)
    gd = groupby(idf,[:Ste,:Cty]) #group by U.S. County FIPS, within the same State FIPS 
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
end

"""
By-airport (`:IATA_Code`) PRE-aggregation, with dumb Euclidean centroid for geographical coordinates.
Input: a FluTE \$name-tracts.dat WITH dsg AP, a la `[:Name,:IATA_Code,:Pop,:LAT,:LNG]`
"""
function aggByAP(idf::DataFrame;make_names=true)
    gd = groupby(idf,[:IATA_Code]) #group by the tract's dsg AP :IATA_Code 
    out = combine(gd, :Pop => sum, :LAT => mean, :LNG => mean, renamecols = false)
    if make_names && "Name" ∉ names(out)
        insertcols!(out, :Name => map( select_mkName(names(out)), eachrow(out)))
    end
    sort!(out,:IATA_Code) #ensure it's sorted by the name/IATA_Code
end

#---AUX::---PARTITION---AND---REVERSE------------------------#
"""Map county names in `pns` to lists of node indices in `ns`."""
function partByCty(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.Ste[ri]==split(x,"~")[1] && 
                        ns.Cty[ri]==split(x,"~")[2],
                        1:nrow(ns)) for x ∈ pns.Name)
end

#consider a unified approach to partBySte and partByAP;
"""Map state names in `pns` to lists of node indices in `ns`."""
function partBySte(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.Ste[ri]==x, 1:nrow(ns)) for x ∈ pns.Ste |> unique)
end

"""Map IATA codes in `pns` to lists of node indices in `ns`."""
function partByAP(ns::DataFrame,pns::DataFrame)
    Dict(x => filter(ri -> ns.IATA_Code[ri]==x, 1:nrow(ns)) for x ∈ pns.IATA_Code |> unique)
end

"""
Reverse-explode a `dict` as returned by `partBy...`, i.e. map indices to their
corresponding parition name.
"""
function revexplPart(dict::Dict)
    ( v => k  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

"""
Reverse-explode a `dict` as returned by `partBy...` and map node names to their
corresponding parition name.
ns MUST have [:Name]
"""
function revexplPart(dict::Dict, ns::DataFrame)
    ( ns.Name[v] => k  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

"""
Reverse-explode a `dict` as returned by `partBy...` and map node names to their
corresponding parition *index*.
ns MUST have [:Name]; ix:: Name => partition_index
"""
function revexplPart(dict::Dict, ns::DataFrame,ix::Dict)
    ( ns.Name[v] => ix[k]  for k ∈ keys(dict) for v ∈ dict[k]) |> Dict
end

end