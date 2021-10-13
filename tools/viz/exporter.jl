"""
Authors: Kara Ignatenko, 2021

Contains functions that help export VegaLite plots to a file.

exporter.jl
2021-09-30 v.0.1: First version
2021-10-13 v.0.2: Remove redundant dispatching by file extension, as VegaLite.jl
                  already has this functionality
"""
module Exporter

using Dates
using VegaLite

export savePlot, figDir

thisPath = splitpath(@__DIR__)
projRoot = thisPath[1:findfirst(isequal("epi-net-m"), thisPath)]
const figDir = joinpath(projRoot..., "fig")

timestamp() = Dates.format(Dates.now(), "yy-mm-dd-HMS")

"""
Export `spec` into a file named `prefix-timestamp.extension`
If `extension` is `"pdf"` or `"csv"`, the spec will be plotted
into the corresponding format. Otherwise a JSON VegaLite spec will be output.
"""
function savePlot(prefix::AbstractString, spec::VegaLite.VLSpec, extension::AbstractString=".html")
    path = joinpath(figDir, prefix * timestamp() * "." * extension)

    save(path, spec)

    path
end


end 