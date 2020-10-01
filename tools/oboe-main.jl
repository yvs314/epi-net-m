#= 
Author: Yaroslav Salii, 2020.

This is Oboe-Mangle, with a view to generate instances for testing
computational methods for networked epidemic models with data from
    (a) FluTE/US 2010 Census Tracts (population, coordinates); github.com/dlchao/FluTE
    (b) US 2019+ Domestic Flights, available from BTS

Reading the FluTE data and transforming it to my input format.
Might separate my input format definitions later.

oboe-main.jl v.0.1: "From scripts to proper code" edition
let's migrate from Jupyter NB, I say.
=#

#reading ingress $name-tract.dat, $name-wf.dat$
#reading ingress; possibly, output too
using CSV
#transforming the data; 
using DataFrames

println("Hello world!")

