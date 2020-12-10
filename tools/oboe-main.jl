"""
Author: Yaroslav Salii, 2020+

Main file for Oboe. Putting all the processing together,
debugging, etc.

Eventually, there'll be a couple debug scenarios and some command-line args parsing

oboe-main.jl
2020-12-10  v.0.0 A Hello, World!

"""

#here's a crutch to load from current path;
#TODO: make a package
push!(LOAD_PATH,pwd())
using Oboe

print(Oboe.callsign)