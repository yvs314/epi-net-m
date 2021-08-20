# epi-net-m
Data Engineering for and Numerical Implementation of Optimal Control of Large-Scale Metapopulation Network SIR

The **Data Engineering** part is found in `/tools` and called “Oboe”; it is written in Julia. Oboe loads the tract populations, commute, and air travel data
from `data/by-tract`; uses this data to generate a travel model aggregated to the given level; and then
outputs a `<name>-trav.dat` matrix containing the total daily travellers between each subdivision 
and a `<name>-init.csv` table containing the initial values for the SIR model.

The **Numerics** is found in `/m-core` and written in Matlab. It uses the _forward-backward sweep_ method to solve a _two-point boundary value problem_ (TPBVP) for a network SIR model instantiated with `<name>-trav.dat` and `<name>-init.csv`. The entry point is `/m-core/sweep.m`. 
A Julia version is 🚧 under construction 🚧

CAVEAT: all data and notebooks are stored with [Git LFS](https://git-lfs.github.com/). If after cloning the repository or downloading its contents, instead of expected file content, you see something like this
```
version https://git-lfs.github.com/spec/v1
oid sha256:9e93547e554054a1678f4863fd62bac1577dd6eea6b2efce0d265b16d6e0f438
size 5208
```
then check your Git LFS installation and try again. 

---
## Usage: Data Engineering

Run the data processing routine from `/tools` as follows:

    julia oboe-cli.jl --fips FIPS... --agg AGG --name NAME [--force]

- `FIPS...` is a space-separated list of [FIPS codes](https://www.nrcs.usda.gov/wps/portal/nrcs/detail/?cid=nrcs143_013696) of the states whose data is to be processed, or `ALL` to process all of the contiguous US, or `NW` to process Washington and Oregon.
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

    # process entire contiguous US, aggregate by state
    julia oboe-cli.jl --fips ALL --agg ste --name USA
---
## Usage: Numerics

Open `/m-core/sweep.m` in MATLAB. Review the _parameters section_ and the requested _instance name_. Run `sweep.m`.

#### Figure export

Semi-automated, look into `./fig` after running the script.
`fAllNodesAbs` tracks per-compartment evolution for all nodes
`fStacked` shows per-node stacked plots of the form `z(t)+s(t)+r(t)`. There's also `m-core/figSimplex.m` that plots each node's trajectory in `(z,s)` simplex.

#### “Time Series” export
Automated, look into `./out` after running the script.
Two sets of `.csv` are produced, giving each simulation _day_'s population of each node's compartments, both in absolute and fractional forms.

#### Compatibility notes

Requires **MATLAB R2019** due to reliance on `tiledlayout` for handling subplots

--- 
## License

The data set and benchmark instances are licensed under GPL v.3.0 since they are in part derived from the [FluTE](https://github.com/dlchao/FluTE) data under the same license. 

🚧 Picking an open license for the code and making this repository public is under construction 🚧
