"""
Author: Kara I, 2021

The CLI entry point for Oboe. This file is responsible for parsing and validating
the CLI arguments and passing them to oboe-main.jl where the actual processing happens.

Run me by typing `include("oboe-main.jl").processOboe(...)` at the julia REPL

oboe-cli.jl
2020-12-10  v.0.0 Initial version
"""

using ArgParse
# using OboeMain
include("oboe-main.jl") # I'm still not entirely sure how packages work in Julia
                        # so using include for now --Kara

#initialize CLI argument parsing
s = ArgParseSettings()

@add_arg_table! s begin
    "--fips"
        nargs = '+'
        help = "FIPS codes for the states that need to be processed, or ALL for all of contiguous US"
        arg_type = String
        required = true
    "--agg"
        nargs = 1
        help = "Aggregation method {tra|cty|ap|ste}"
        arg_type = String
        required = true
    "--name"
        nargs = 1
        help = "The name of the dataset to be used, must not contain numbers"
        arg_type = String
        required = true
end

function validate_args(agg, name)
    if agg ∉ ["tra", "cty", "ap", "ste"]
        throw(DomainError(agg, "The argument for agg must be {tra|cty|ap|ste}"))
    end
    if any(ch -> '0' <= ch <= '9', name)
        throw(DomainError(name, "The name must not contain numbers"))
    end
end

useall = false # true if all contiguous US states should be used
useNW = false # true if the a~NW-tracts.dat dataset should be used instead of usa-tracts.dat

#Parse and process args
parsed_args = parse_args(ARGS, s)
#Get the target FIPS codes from the args, or use all states is ALL is given instead
fips = parsed_args["fips"]
if fips == ["ALL"]
    useall = true
elseif fips == ["NW"]
    useNW = true
    fips = ["41", "53"]
end

#Get the aggregation mode and the filename
agg = parsed_args["agg"][1]
name = uppercase(parsed_args["name"][1])

validate_args(agg, name)

#perform the data processing routine using the given args

if useall
    OboeMain.processOboe(name, agg)
elseif useNW
    OboeMain.processOboe(name, agg, fips=fips, useNW=true)
else
    OboeMain.processOboe(name, agg, fips=fips)
end
