"""
Authors:  Yaroslav Salii, 2021
          Kara Ignatenko, 2021

Contains VegaLite spec-generating functions to be used in visualizations.

vega-specs.jl
2021-09-22 v.0.1: First modular version
2021-09-28 v.0.2: Add more specs for network visualization
2021-10-09 v.0.3: Add specs for the average infected rate line plot
2021-10-13 v.0.4: Drop the leading 0's of FIPS codes to resolve the inconsistency
                  between how FIPS codes are represented in Oboe vs. VegaDatasets
2021-10-24 v.0.5: Add a spec for control effort choropleth
"""
module VegaSpecs

using VegaLite, VegaDatasets, DataFrames

us10m = dataset("us-10m")

"""Utility function to drop the leading zero of a FIPS code."""
function dropLeadingZero(str::String)
    if startswith(str, "0")
        str[2:length(str)]
    else
        str
    end
end

normalize(arr::AbstractArray) = map(dropLeadingZero, arr)

"""
Plot side-by-side chloropleths comparing the Z and Z0 values in the solution `sol`
for a given day `day`, with `median` as the centre of the colour scale.
"""
function pltZOptVsNullByCty(sol::DataFrame, day::Int, median::Float64)
    # helper functions to get the per names of the Z and Z0 columns for a given day
    day2col(day::Int)  = "Z" * string(day)
    day2col0(day::Int) = day2col(day) * "_NULL"

    @vlplot(
	repeat={column=[day2col(day), day2col(day)*"_NULL"]},
    ) + 
    @vlplot(
        :geoshape,
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
                data=select(sol, ["id", day2col(day), day2col0(day)]), #INSERT DF HERE
                key=:id,
                fields=[day2col(day), day2col0(day)]
            }
        }],
        projection={
            type=:albersUsa
        },
        encoding={
            color={
                field={repeat=:column},
                type=:quantitative,
                scale={
                    domainMid= median, #median no. infected
                    scheme=:reds
                }
            }
        }
    )
end

"""
Plot a control effort by-county choropleth map for the given `day`
The control effort u has a blue-to-purple scheme and [0,1] scale
"""
function pltCtrlByCty(sol::DataFrame, day::Int)
    day2col_u(day::Int) = "u" * string(day) #picking the control column name for the day
    u_today = day2col_u(day) #compute the requisite column name
    @vlplot(
        :geoshape,
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
                data=select(sol, ["id", u_today]), #INSERT DF HERE
                key=:id,
                fields=[u_today]
            }
        }],
        projection={
            type=:albersUsa
        },
        encoding={
            color={
                field= u_today,
                type=:quantitative,
                scale={
                    domainMin = 0,
                    domainMax = 1, #median no. infected
                    scheme=:bluepurple
                }
            }
        }
    )
end

# """
# ISSUES: botches the color scale properties
# Plot 3 side-by-side chloropleths comparing the u, Z, and Z_Null values in the solution `sol`
# for a given day `day`.
# For Z and Z_Null, the `median` as the centre of the Z and Z_Null colour scale.
# The control effort u is displayed to the left, with green scheme and [0,1] scale
# """
# function pltCtrl_ZOpt_ZNull_ByCty(sol::DataFrame, day::Int, median::Float64)
#     @vlplot() + [pltCtrlByCty(sol,day), pltZOptVsNullByCty(sol,day,median)] 
# end

"""
Plot the average infection rate Z (for both the optimal case
and the null case) and the control effort U as functions
of time.
`long_avgc` must have `[:z_avg, :zNull_avg, :u_avg]`
"""
pltAvgInfdCtrl(long_avgc::DataFrame) = @vlplot(
    data = long_avgc,
    :line,
    x=:day,
    y=:value,
    color={
        :symbol, 
        scale={
            domain=["zNull_avg","z_avg","u_avg"],
            range = ["#E31A1C","#FB9A99","#33A02C"] #red, pink, green
            #scheme=:paired
        },
        #legend = nothing
        legend = {
            orient = "top-right",
            #values=["lorem","ipsum","sit"]
            title = "variable",
        }
    },
)


pltCanvasSte(fips::Vector{String}) = @vlplot(
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
            filter={field =:id,oneOf = normalize(fips)}
            }],
    projection={type=:albersUsa},
)

pltCanvasCty(ctys::DataFrame) = @vlplot(
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

pltBoundarySte(fips::Vector{String}) = @vlplot(
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
            filter={field =:id,oneOf = normalize(fips)}
            }],
    projection={type=:albersUsa},
)

pltDotsTracts(ns::DataFrame) = @vlplot(
    :circle,
    data=ns,
    projection={type=:albersUsa},
    longitude="LNG:q",
    latitude="LAT:q",
    tooltip="Name:n",
    size={value=2},
    color={value=:black})

pltDotsDsgAPs(dsgAPs::DataFrame) = @vlplot(
    :circle,
    data=dsgAPs,
    projection={type=:albersUsa},
    longitude="LNG:q",
    latitude="LAT:q",
    tooltip="IATA_Code:n",
    size={value=60},
    color={value=:steelblue}
)

function pltTractsAndAPsBySte(ns::DataFrame, dsgAPs::DataFrame, fips::Vector{String})
    @vlplot() + pltCanvasSte(fips) + pltDotsTracts(ns) + pltDotsDsgAPs(dsgAPs)
end

function pltTractsAndAPsByCty(ns::DataFrame, dsgAPs::DataFrame, fips::Vector{String}, ctys::DataFrame; borders=true::Bool)
    @vlplot() + pltCanvasCty(ctys) +  pltDotsTracts(ns) + pltDotsDsgAPs(dsgAPs) + pltBoundarySte(fips) 
end

pltPopCty(ctys::DataFrame) = @vlplot(
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

function pltPopCtyWithBorders(ctys::DataFrame, fips::Vector{String})
     @vlplot() + pltPopCty(ctys) + pltBoundarySte(fips)
end

end
