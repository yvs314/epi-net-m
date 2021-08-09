"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

This is Oboe-Mangle, or maybe Arshyn, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS
"""
module Oboe
using FromFile



const callsign = "This is Oboe v.1.0" * begin
    branch = readchomp(`sh -c "git branch --show-current 2> /dev/null || echo master"`)
    if branch == "master"
        ""
    else
        "-" * branch
    end
end

@from "io.jl"          using IO
@from "airports.jl"    using Airports
@from "tracts.jl"      using Tracts
@from "aggregation.jl" using Aggregation
@from "travel.jl"      using Travel

end