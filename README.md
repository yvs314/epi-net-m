# epi-net-m
Simulation Tools for Networked SIR+

The Oboe data processing routine will load the tract populations, commute and air travel data
from `data/by-tract`, use this data to generate a travel model aggregated to the given level and
output a `.dat` matrix containing the total daily travellers between each subdivision 
and a `.csv` table containing the initial values for the SIR model.

## Usage: Data Processing

Run the data processing routine as follows:

    julia oboe-cli.jl --fips FIPS --agg AGG --name NAME [--force]

- `FIPS` is a space-separated list of [FIPS codes](https://www.nrcs.usda.gov/wps/portal/nrcs/detail/?cid=nrcs143_013696) of the states whose data is to be processed, or `ALL` to process all of the contiguous US, or `NW` to process Washington and Oregon.
- `AGG` is the desired aggregation level, which can be one of:
    - `tra` for census tract
    - `cty` for county
    - `ap` for airport
    - `ste` for state
- `NAME` is the prefix for output files. Will be treated as all-uppercase. Must not contain dashes (`-`) or underscores (`_`).
- Include `--force` to allow overwriting existing output files with the same prefix.

### Examples

    # process Washington and Oregon, aggregate to tract level
    julia oboe-cli.jl --fips NW --agg tra --name NW

    # process California, aggregate by airport, force overwrite
    julia oboe-cli.jl --fips 06 --agg ap --name CA --force

    # process NY, NJ and CT, aggregate by county, force overwrite
    julia oboe-cli.jl --fips 09 34 36 --agg cty --name TRI --force

## Usage

1. Download/clone the repository

2. Load main_script.m into MATLAB and run it.

3. Admire the produced figures.

4. Stand in awe of the produced tables.

## Figure export

Automated, look into `./out` after running the script.
`fAllNodesAbs` tracks per-compartment evolution for all nodes
`fStacked` shows per-node stacked plots of the form `i(t)+s(t)+r(t)`

## “Time Series” export
Automated, look into `./out` after running the script.
Two `.csv` are produced, giving each simulation _day_'s population of each node's compartments, both in absolute and fractional forms.

## Compatibility notes

Requires **MATLAB R2019** due to reliance on `tiledlayout` for handling subplots

## License

The data set/benchmark instances are licensed under GPL v.3.0 since they are in part derived from the FluTE data under the same license.
