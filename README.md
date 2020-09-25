# epi-net-m
Simulation Tools for Networked SIR+

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
