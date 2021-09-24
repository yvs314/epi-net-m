"""
Authors: Kara Ignatenko, 2021

Contains VegaLite spec-generating functions to be used in PLuto.jl visualizations.

vega-specs.jl
2021-09-22 v.0.1: First modular version
"""
module VizSpecs

using VegaLite, VegaDatasets, DataFrames

us10m = dataset("us-10m")

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

end