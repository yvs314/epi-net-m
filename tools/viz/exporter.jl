"""
Authors: Kara Ignatenko, 2021

Contains functions that help export VegaLite plots to a file.

exporter.jl
2021-09-30 v.0.1: First version
"""
module Exporter

using Dates
using VegaLite

export savePlot

thisPath = splitpath(@__DIR__)
projRoot = thisPath[1:findfirst(isequal("epi-net-m"), thisPath)]
const figDir = joinpath(projRoot..., "fig")
const MIMEtypes = Dict(
    "pdf" => MIME"application/pdf",
    "svg" => MIME"image/svg+xml"
)
const defaultType = MIME"application/vnd.vegalite.v4+json"

timestamp() = Dates.format(Dates.now(), "yy-mm-dd-HMS")

"""
Export `spec` into a file named `prefix-timestamp.extension`
If `extension` is `"pdf"` or `"csv"`, the spec will be plotted
into the corresponding format. Otherwise a JSON Vega spec will be output.
"""
function savePlot(prefix::AbstractString, spec::VegaLite.VLSpec, extension::AbstractString)
    path = joinpath(figDir, prefix * timestamp() * "." * extension)

    try
        io = open(path, "w")
        show(io, get(MIMEtypes, extension, defaultType), spec)
        close(io)
    catch e
        "Error while writing: " * string(e)
    end

    path
end


end 