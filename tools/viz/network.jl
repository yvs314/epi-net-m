module Network

using DataFrames, Mux, FromFile
using VegaLite: VLSpec
# using HTTP: Request, Response
@from "../Oboe/Oboe.jl" using Oboe
@from "./vega-specs.jl" using VegaSpecs


# stateCache = Dict{String, State}()


### Before the server even starts, we load and preprocess everything

nsRaw = Oboe.rdWholeUS()
APs = Oboe.getProcessedAPs(Oboe.rdBTS() |> Oboe.grpBTS, Oboe.rdAPs())
ns = Oboe.assignDsgAPs(nsRaw,APs)
println("Network information loaded")

"""
Return a NamedTuple where `tracts` is all the nodes in the selected states,
`APs` is the list of designatedd airports in those states and `ctys` is all
the counties
"""
function states(fips::Vector)
    normFIPS = map(Oboe.normalizeOneDigitFIPS, fips)
    tracts = Oboe.censorFluteTractByFIPS(ns, normFIPS)

    dsgAPcodes = ns.IATA_Code |> unique
    dsgAPs = filter(r -> r.IATA_Code ∈ dsgAPcodes,eachrow(APs)) |> DataFrame

    # cheapest way of filtering for _relevant_ county IDs in plots
    ctys = Oboe.aggByCty(tracts)
    insertcols!(ctys,1,:id => 
       map(r-> parse(Int,r.Ste*"000") + parse(Int,r.Cty),eachrow(ctys)) )

    return (tracts=tracts, APs=dsgAPs, ctys=ctys, fips=normFIPS)
end

stringify(substrs::Vector{SubString{String}}) = map(String, substrs)
parseFIPS(str::AbstractString) = map(Oboe.normalizeOneDigitFIPS, split(str, ",")) |> stringify
getStates(req::Dict) = states(parseFIPS(req[:params][:fips]))

### API functions below ###

netHandler = mux(
    route("/bySte/:fips", req -> begin
        stes = getStates(req)
        VegaSpecs.pltTractsAndAPsBySte(stes.ns, stes.dsgAPs, stes.fips)
    end),
    route("/byCty/:fips", req -> begin
        stes = getStates(req)
        VegaSpecs.pltTractsAndAPsByCty(stes.ns, stes.dsgAPs, stes.fips, stes.ctys)
    end),
    route("/pop/:fips", req -> begin
        stes = getStates(req)
        VegaSpecs.pltPopCtyWithBorders(stes.ctys, stes.fips)
    end),
    Mux.notfound()
)

export netHandler

end