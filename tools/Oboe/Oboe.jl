"""
Authors: Yaroslav Salii, 2020+
         Kara Ignatenko, 2021

This is Oboe-Mangle, or maybe Arshyn, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS

Reading the FluTE data and transforming it to my input format.
Might separate my input format definitions later.

Oboe.jl v.0.1: "From scripts to proper code" edition
        v.0.2: switch to `module`
        v.0.3: got as far as reading FluTE tracts
        v.0.4: added basic by-county and by-state IV aggregation
        v.0.5: join (intersect) BTS with OpenFlights
        v.0.6: added designated AP assignment
'20-12-18   v.0.7: wrote pairs-to-matrix xform by hand (vs. `unstack`, which was unpredictable)
'20-12-30   v.0.8: added a beta node-to-node daily air passenger computation
'21-01-03   v.0.8.1: a half-baked Main(), look in bit bucket. Tested on 3K by-county!
'21-01-22   v.0.9: air travel and commute matrices are in, tested on 2K NW (by-tract)
'21-02-02   v.0.9.1: hotfix, removing deprecated `by`, `names!`, and the like (half-done)
'21-02-02   v.0.9.2: removed all the deprecated stuff, and it all works again
'21-02-14   v.0.9.3: add node partitioning, streamline mkFlightMx2() and mkPsgMx2() with aux
'21-02-14   v.0.9.4: add partition-to-partition flight matrix
'21-02-18   v.0.9.5: add partition-to-partition commute matrix
'21-03-31   v.0.9.6: min,median,max for nonzero elements in talkDnsy
            v.0.A: add by-AP agg with recomputation bypass in mkPsgMx()
'21-07-14   v.0.A.1: add checkIfFilesExist and censorFluteTractByFIPS
'21-07-19   v.1.0: fixed issue where unmatched tract IDs led to a crash
'21-08-10   v.1.1: significantly improved performance of mkPsgMx
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