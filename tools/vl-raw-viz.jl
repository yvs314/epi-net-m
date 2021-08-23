"""
A VegaLite.jl-based plotting routine, 
written for line-by-line execution through Julia-in-VS-Code
2021-08-23 v.0.0 can draw APs, states, and census tracts
"""

using VegaLite,VegaDatasets
using DataFrames, CSV, Statistics, Revise

pwd()

include("./Oboe/Oboe.jl")

rawBTS = Oboe.rdBTS()

#println(Oboe.callsign)

fips=["41","53"]

nsRaw = Oboe.censorFluteTractByFIPS(Oboe.rdWholeUS(), fips)

APs = Oboe.getProcessedAPs(Oboe.rdBTS() |> Oboe.grpBTS, Oboe.rdAPs())

ns = Oboe.assignDsgAPs(nsRaw,APs)

#byctyRaw = Oboe.aggByCty(ns)


dsgAPcodes = ns.IATA_Code |> unique

dsgAPs = filter(r -> r.IATA_Code ∈ dsgAPcodes,eachrow(APs)) |> DataFrame

#====PLOT==TRACTS==APs==COUNTIES===========#

us10m = dataset("us-10m")

pltCanvasSte = @vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:white
    },
    data={
        values=us10m,
        format={
            type=:topojson,
            feature=:states
        }
    },
    transform=[{
            filter={field =:id,oneOf = fips}
            }],
    projection={type=:albersUsa},
)

# pltCanvasCty = @vlplot(
#     mark={
#         :geoshape,
#         fill=:lightgray,
#         stroke=:white
#     },
#     data={
#         values=us10m,
#         format={
#             type=:topojson,
#             feature=:counties
#         }
#     },
#     transform=[{
#             filter={field =:id,oneOf = fips}
#             }],
#     projection={type=:albersUsa},
# )

pltDotsTracts = @vlplot(
    :circle,
    data=ns,
    projection={type=:albersUsa},
    longitude="LNG:q",
    latitude="LAT:q",
    tooltip="Name:n",
    size={value=2},
    color={value=:black})

pltDotsDsgAPs = @vlplot(
    :circle,
    data=dsgAPs,
    projection={type=:albersUsa},
    longitude="LNG:q",
    latitude="LAT:q",
    tooltip="IATA_Code:n",
    size={value=60},
    color={value=:steelblue}
)

#NB! the first empty `@vlplot()` is necessary to use the concat + syntax
pltM = @vlplot() + pltCanvasSte + pltDotsTracts + pltDotsDsgAPs

save("../fig/test.pdf",pltM)

#pltM2 = @vlplot() + pltCanvasCty + pltDotsTracts + pltDotsDsgAPs