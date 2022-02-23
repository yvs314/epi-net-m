# epi-net-m

## A Benchmarking Framework for Optimal Control over Network Dynamic Systems based on a Metapopulation Epidemic Model
With **Data Engineering** to generate a **Benchmark**, a proof-of-concept **Numerical Solver** for an optimal control problem over Network SIR, and **Visualization Tools** to map infection and control effort on [choropleths](https://en.wikipedia.org/wiki/Choropleth_map). The _optimal control problem_ is finite-horizon Bolza over a nonlinear system with per-node _social distancing_ controls and time-discounted _running cost_, square in control and infections. It generalizes [(El Ouardighi, Khmeltinsky, Sethi, 2022)](https://doi.org/10.1111/poms.13641) to the network dynamic system case, with omission of health infrastructure tracking.



The **Data Engineering** part is in `/tools` and called “Oboe”; it is written in Julia. Oboe loads the tract populations, commute, and air travel data
from `/data/by-tract`; uses this data to generate a travel model aggregated to the given level; and then
outputs a `<name>-trav.dat` matrix containing the total daily travellers between each subdivision 
and a `<name>-init.csv` table containing the “patient zero” initial values for the SIR model. These serve as a benchmark set.

The **Benchmark** is in `/data/by-tract`, with some snapshots in [Releases](https://github.com/yvs314/epi-net-m/releases); these have between 2 and 9110 nodes. An instance `<name>` is made of two files,
- `<name>-init.csv` with one row for each of the `n` nodes listing the populations and initial values of _susceptible_, _infected_, and _recovered_,  
- `<name>-trav.dat` with `n×n` matrix of daily travelers between the nodes

The **Numerics** is in `/m-core` and written in Matlab. It uses the _forward-backward sweep_ method to solve a _two-point boundary value problem_ (TPBVP) for a network SIR model instantiated with `<name>-trav.dat` and `<name>-init.csv`. The entry point is `/m-core/sweep.m`. 
🚧 A Julia version may or may not be under construction 🚧

The **Visualization** is in `/tools/viz` and based on [VegaLite](https://vega.github.io/vega-lite/)'s Julia interface. The VegaLite schemas generate _county-level_ choropleths for infection incidence and control effort.


CAVEAT: all data and Jupyter notebooks are stored with [Git LFS](https://git-lfs.github.com/). If after cloning the repository or downloading its contents, instead of expected file content, you see something like this
```
version https://git-lfs.github.com/spec/v1
oid sha256:9e93547e554054a1678f4863fd62bac1577dd6eea6b2efce0d265b16d6e0f438
size 5208
```
then your Git LFS installation did not work. Get the benchmark from [Releases](https://github.com/yvs314/epi-net-m/releases) if you are not in the mood for Git LFS.

---
## Citation
If you use this software, please cite this repository and [(Salii, 2022)](https://doi.org/10.1007/978-3-030-93413-2_17
)
```
@inproceedings{salii2022benchmarking,
	author={Salii, Yaroslav V.},
	editor={Benito, Rosa Maria and Cherifi, Chantal and Cherifi, Hocine and Moro, Esteban and Rocha, Luis M. and Sales-Pardo, Marta},
	title={Benchmarking Optimal Control for Network Dynamic Systems with Plausible Epidemic Models},
	booktitle={Complex Networks {\&} Their Applications {X}},
	year={2022},
	publisher={Springer International Publishing},
	address={Cham},
	pages={194--206},
	isbn={978-3-030-93413-2}
}
```

---
## Usage: Data Engineering

Run the data processing routine from `/tools` as follows:

    julia oboe-cli.jl --fips FIPS... --agg AGG --name NAME [--force]

- `FIPS...` is a space-separated list of [FIPS codes](https://www.nrcs.usda.gov/wps/portal/nrcs/detail/?cid=nrcs143_013696) of the states whose data is to be processed, or `ALL` to process all of the contiguous US. 
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

### Julia API

To access the Oboe API from within Julia, execute the following:

```julia
using FromFile
@from "(path to)/Oboe/Oboe.jl" import Oboe
```
- If running from a Julia file or Jupyter notebook, the path should be relative to the location of said file.
- If running from a REPL, the path should be relative to the working directory of the REPL. In particular, when using VS Code's "execute in REPL" functionality, the working directory is set to the project root, so the path should be `./tools/Oboe/Oboe.jl`

This should bring every API function into the `Oboe` namespace, e.g. `Oboe.mkPsgMx`.

---
## Usage: Numerics

Open `/m-core/sweep.m` in MATLAB. Review the _parameters section_ and the requested _instance name_. Run `sweep.m`.

### Solution export
Automated, look into `./out` after running `sweep.m`.
Four sets of `.csv` are produced, giving each simulation _day_'s population of each node's compartments, both in absolute and fractional forms for _optimal control_ and _null control_. The _per-node_ **control effort** is exported in `<name>-frac.csv` as `uX` columns, where `X` is the simulation day's number. In addition, a `<name>-log.csv` is emitted, which describes the forward-backward sweep iterations.


### Figures 
- `figStacked.m` _stacked plot_ of `z+s+r` for a given node
- `figSimplex.m` all the nodes' trajectories in the `(z,s)` simplex (_% infected_, _%susceptible_)
- `figTrajectory.m` all the nodes' trajectories for a given _compartment_: susceptible `s`, infected `z`, or recovered `r`

--- 
## Usage: Visualization
See the VegaLite-based routines in [`/tools/viz`](https://github.com/yvs314/epi-net-m/tree/master/tools/viz). 

There are two Pluto.jl-based notebooks, which work for county-level aggregation and provide `.svg` export:
- `network-explorer.jl` Displays the _designated airports_ as computed with Oboe, the census tracts, and also per-county populations
- `solution-explorer.jl` From `.csv` _solution files_ displays 
    - the absolute numbers of infected per county, side-by-side in **optimal control** vs **null** control
    - the **control effort** per-county
    - average **control effort** plot (average over populations)

The VegaLite **schemas** used in the above two are in `vega-specs.jl` and can be used independently, e.g. through Jupyter Notebook or Julia-for-VS-Code.

---
## Acknowledgements

The data set and benchmark instances are in part derived from the [FluTE](https://github.com/dlchao/FluTE) data, coupled with U.S. domestic carrier air travel data from [U.S. Bureau of Transportation Statistics](https://www.transtats.bts.gov/) and airport information from [Openflights repository](https://github.com/jpatokal/openflights). 

**Yaroslav Salii** @yvs314 is the principal author, who designed the original version of the Oboe data processing routine and the MATLAB implementation of the [Forward-Backward Sweep](https://github.com/yvs314/epi-net-m/blob/8584d09125a2250032ff8300365daa92fe3941e4/m-core/sweep.m) numerical solution method for Network Metapopulation SIR Epidemic Model with Social-Distancing Optimal Control, and VegaLite-based visualizations.

**Kara Ignatenko** @karaign implemented the Oboe command-line interface, significantly improved the performance of air travel network generator [mkPsgMx](https://github.com/yvs314/epi-net-m/blob/8584d09125a2250032ff8300365daa92fe3941e4/tools/Oboe/travel.jl), implemented the [FromFile.jl](https://github.com/Roger-luo/FromFile.jl)-based modular version of the Oboe data processing routine, and Pluto.jl-based _solution explorers_.

**Rinel Foguen Tchuendom** and **Shuang Gao** @sigmagao jointly contributed an early version of the [Euler method solver](https://github.com/yvs314/epi-net-m/blob/8584d09125a2250032ff8300365daa92fe3941e4/m-core/old-eulerkron.m) for Network Metapopulation SIR Epidemic Model in Kronecker product notation.

This software was written when the authors were with _McGill University_.