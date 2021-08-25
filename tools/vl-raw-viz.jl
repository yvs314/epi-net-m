"""
A VegaLite.jl-based plotting routine, 
written for line-by-line execution through Julia-in-VS-Code

This one is for _a priori_ figures: 
(a) census tracts with designated APs and state boundaries (b) county populations


2021-08-23 v.0.0 can draw APs, states, and census tracts, pltM2
2021-08-23 v.0.1 can also draw county population, pltP
"""

using VegaLite,VegaDatasets
using FromFile, DataFrames, CSV, Statistics, Revise

@from "./Oboe/Oboe.jl" import Oboe
#include("./Oboe/Oboe.jl")



#====PREP==THE==DATA=================#

fips = ["41","53"]

nsRaw = Oboe.censorFluteTractByFIPS(Oboe.rdWholeUS(), fips)

APs = Oboe.getProcessedAPs(Oboe.rdBTS() |> Oboe.grpBTS, Oboe.rdAPs())

ns = Oboe.assignDsgAPs(nsRaw,APs)

#cheapest way of filtering for _relevant_ county IDs in plots
ctys = Oboe.aggByCty(ns)
insertcols!(ctys,1,:id => 
    map(r-> parse(Int,r.Ste*"000") + parse(Int,r.Cty),eachrow(ctys)) )

dsgAPcodes = ns.IATA_Code |> unique

dsgAPs = filter(r -> r.IATA_Code ∈ dsgAPcodes,eachrow(APs)) |> DataFrame

#====PLOT==TRACTS==APs==COUNTIES===========#

#load the US 1:10^6 GeoJSON shapes

us10m = dataset("us-10m")

pltCanvasSte = @vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:white,
        strokeWidth=1
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

pltCanvasCty = @vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:white
    },
    data={
        values=us10m,
        format={
            type=:topojson,
            feature=:counties
        }
    },
    transform=[{
            filter={field =:id,oneOf = ctys.id}
            }],
    projection={type=:albersUsa},
)

pltBoundarySte = @vlplot(
    mark={
        :geoshape,
        stroke=:white,
        filled = false,
        strokeWidth=3,
        strokeOpacity=0.2
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
pltM2 = @vlplot() + pltCanvasCty +  pltDotsTracts + pltDotsDsgAPs + pltBoundarySte 

print(pwd())
save("../fig/NWtracts_aps.pdf",pltM)
save("../fig/NWtracts_aps2.pdf",pltM2)


#====COUNTY===POP===PLOT=======#

pltPopCty = @vlplot(
    :geoshape,
   # width=500, height=300,
    data={
        values=us10m,
        format={
            type=:topojson,
            feature=:counties
        }
    },
    transform=[{
        lookup=:id,
        from={
            data=ctys, #INSERT DF HERE
            key=:id,
            fields=["Pop"]
        }
    }],
    projection={
        type=:albersUsa
    },
    encoding={
        color={
            field="Pop",
            type=:quantitative,
            # type=:log,
            scale={
                domainMid= 44479, #median county pop
                scheme=:blues
                }
        }
    }
    
)

pltP = @vlplot() + pltPopCty + pltBoundarySte

pwd() |> print
save("../fig/NWctyPop2.pdf",pltP)

