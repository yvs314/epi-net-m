# epi-net-m
Simulation Tools for Networked SIR+

## Usage

1. Download/clone the repository

2. Load main_script.m into MATLAB and run it.

3. Admire the produced figures.

## Figure export

At this time, only manual. Suggested location: `./fig`, not tracked by the version control.
`fAllNodesAbs` tracks per-compartment evolution for all nodes
`fStacked` shows per-node stacked plots of the form `i(t)+s(t)+r(t)`

Use `print`, `exportgraphics`, or something else on these

## Compatibility notes

Requires **MATLAB R2019** due to reliance on `tiledlayout` for handling subplots
