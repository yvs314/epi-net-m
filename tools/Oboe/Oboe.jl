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

const callsign="This is Oboe v.1.0"

@from "io.jl"          import IO
@from "airports.jl"    import Airports
@from "tracts.jl"      import Tracts
@from "aggregation.jl" import Aggregation
@from "travel.jl"      import Travel

end