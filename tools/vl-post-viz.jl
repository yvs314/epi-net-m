"""
A VegaLite.jl-based plotting routine, 
written for line-by-line execution through Julia-in-VS-Code

This one is mostly for _a posteriori_ figures: 
(a) initial conditions/patient zero location
(b) infected at last day with null control
(c) infected at last day with optimal control

Yaroslav Salii, 2021
2021-08-25 v.0.0 all plots automated
2021-10-04 v.0.1 adding avg infected rate and control line plot
"""

using VegaLite,VegaDatasets
using FromFile, DataFrames, CSV, Statistics, Revise


thisPath = splitpath(@__DIR__)
projRoot = thisPath[1:findfirst(isequal("epi-net-m"), thisPath)]

const figDir = joinpath(projRoot...,"fig") #where to put the resulting plots
const iDir = joinpath(projRoot...,"out") #here be the simulator's outputs

#=
note: @__DIR__'s output is valid (via VS Code), pwd()'s not, as far as @from is concerned
will track this issue in `save` functions, hopefully circumvented by absolute path in `figDir`
=#
@__DIR__
pwd()

#const oboePath = joinpath(projRoot...,"tools","Oboe","Oboe.jl")
#@from "./Oboe/Oboe.jl" import Oboe

#====PREP==THE==DATA=================#
slnName="NWcty_75"
#paths to solutions in fractions, optimal control and null-control
#slnOptPath=joinpath(iDir,slnName*"-frac.csv")
#slnNullPath=joinpath(iDir,slnName*"-frac0.csv")

"""
Read the solutions into dataframes and stuff them into a named tuple 
:f,:f0 are *fractional* per-node per-day s/z/r
:a,:a0 are *absolute* per-node per-day S/Z/R 
:avgc are per-day average infected fraction z, zNull, and average optimal control effort u
"""
function rdSolutions(;iName = slnName, slnDir = iDir )
    slnSuffs = ["-frac.csv","-frac0.csv","-abs.csv","-abs0.csv","-avg.csv"]
    slns= map( p -> CSV.read(p,DataFrame), (joinpath(iDir,slnName * suff) for suff in slnSuffs))
    return NamedTuple([:f,:f0,:a,:a0,:avgc] .=> slns) 
end

#do read all the solutions
ss = rdSolutions()

#=
"-avg.csv" has 3 columns: z_avg; zNull_avg; u_avg, and per-day rows. Gotta fix that in MATLAB. 
also, I should isolate the fundamental plot specs as functions that take the data ("INSERT DF HERE")
=#

a0Median = ss.a0.Z180 |> median
a0Max = ss.a0.Z180 |> maximum

#====PLOT==ABS==INFECTED=======#

#load the US 1:10^6 GeoJSON shapes

us10m = dataset("us-10m")

#plot the no-control case
pltZNull180 = @vlplot(
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
            data=ss.a0, #INSERT DF HERE
            key=:id,
            fields=["Z180"]
        }
    }],
    projection={
        type=:albersUsa
    },
    encoding={
        color={
            field="Z180",
            type=:quantitative,
            scale={
                domainMid= a0Median, #median no. infected
                scheme=:reds
                }
        }
    }
)

#plot the optimal control case
pltZOpt180 = @vlplot(
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
            data=ss.a, #INSERT DF HERE
            key=:id,
            fields=["Z180"]
        }
    }],
    projection={
        type=:albersUsa
    },
    encoding={
        color={
            field="Z180",
            type=:quantitative,
            scale={
                domainMid= a0Median, #median no. infected
                domainMax = a0Max,
                scheme=:reds
                }
        }
    }
)

#plot the initial values
pltIVs_ = @vlplot(
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
            data=ss.a, #INSERT DF HERE
            key=:id,
            fields=["I_i"]
        }
    }],
    projection={
        type=:albersUsa
    },
    encoding={
        color={
            field="I_i",
            type=:quantitative,
            scale={
                scheme=:reds,
                }, 
            legend = nothing, #“nothing” in VegaLite.jl, instead of Vega-Lite's “null”
        }
    }
)

#note: reused from vl-pre-viz.jl, candidate for isolation
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
            filter={field =:id,oneOf = fips} #insert state(?) FIPS
            }],
    projection={type=:albersUsa},
)


#NB! the first empty `@vlplot()` is necessary to use the concat + syntax
pltOpt = @vlplot() + pltZOpt180 + pltBoundarySte 
pltNull = @vlplot() + pltZNull180 + pltBoundarySte 
pltIVs = @vlplot() + pltIVs_ + pltBoundarySte 

plotPaths=[joinpath(figDir,slnName * figSuff *".pdf") 
    for figSuff ∈ ["-Z180opt","-Z180null","-IVs"]]

save(plotPaths[1],pltOpt)
save(plotPaths[2],pltNull)
save(plotPaths[3],pltIVs)
#save all the figures
#map(p -> save(p[1],p[2]),(plotPaths .=> [pltOpt,pltNull,pltIVs]))

#=========BIT====BUCKET===============#
#note: reused from vl-pre-viz.jl, candidate for isolation
# pltCanvasSte = @vlplot(
#     mark={
#         :geoshape,
#         fill=:lightgray,
#         stroke=:white,
#         strokeWidth=1
#     },
#     data={
#         values=us10m,
#         format={
#             type=:topojson,
#             feature=:states
#         }
#     },
#     transform=[{
#             filter={field =:id,oneOf = fips}
#             }],
#     projection={type=:albersUsa},
# )
# #note: reused from vl-pre-viz.jl, candidate for isolation
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
#             filter={field =:id,oneOf = ctys.id}
#             }],
#     projection={type=:albersUsa},
# )