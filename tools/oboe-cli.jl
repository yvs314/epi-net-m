"""
Author: Kara I, 2021

The CLI entry point for Oboe. This file is responsible for parsing and validating
the CLI arguments and passing them to oboe-main.jl where the actual processing happens.

Run me as follows:
julia oboe-cli.jl --fips FIPS --agg AGG --name NAME [--force]

oboe-cli.jl
2021-07-10  v.0.0 Initial version
2021-07-14  v.0.1 Added support for --force
2021-07-15  v.0.1.1 YS: clarify help and validation error text
"""

using FromFile 
using ArgParse
@from "oboe-main.jl" import OboeMain

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
        help = "Pick aggregation type {tra|cty|ap|ste}"
        arg_type = String
        required = true
    "--name"
        nargs = 1
        help = "Output file names prefix, UPPERCASE. No hyphens - or underscores _"
        arg_type = String
        required = true
    "--force"
        action = :store_true
        help = "Overwrite existing output files"
end

function validate_args(agg, name)
    if agg ∉ ["tra", "cty", "ap", "ste"]
        throw(DomainError(agg, "Unknown aggregation type. Choose {tra|cty|ap|ste}"))
    end
    if any(char -> char ∈ ['-','_'], name)
        throw(DomainError(name, "Found a hyphen - or underscore _ in output prefix. These are not allowed"))
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
force = parsed_args["force"]

validate_args(agg, name)

#perform the data processing routine using the given args

if useall
    OboeMain.processOboe(name, agg, force=force)
else
    OboeMain.processOboe(name, agg, fips=fips, useNW=useNW, force=force)
end
