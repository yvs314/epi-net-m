"""
This module is responsible for handling input and output, i.e.
reading tract and travel data from disk and outputting travel matrices.
"""
module IO
using DataFrames
#CSV: IO FluTE & my, DelimitedFiles: writing the matrices (.dat)
using CSV,DelimitedFiles

#====BASE===FILENAMES==TYPES==DATA=STRUCTURES=====#
#all you need to know about input and output of instances
struct NamingSpec #all fields are String, don't say I didn't warn you
    # e.g. sep="-", sSep="_"; $name_$size-$suff
    sep::String #to the right of last $sep is filename suffix, to the left is the instance name
    sSep::String
    #read from ifDir, write to ofDir
    ifDir::String
 #   ifDirAir::String
    ofDir::String
    # what's after instance's name in its Initial Values file name
    fltInitSuff::String
    myInitSuff::String
end

#the naming conventions I am going to use
#as well as input and output directories
global const fn=NamingSpec("-","_"
    ,joinpath("..","data","by-tract","flute")
    ,joinpath("..","data","by-tract")
    ,"tracts.dat","init.csv")

function checkIfFilesExist(iname::String; ofdir=fn.ofDir)
    dir = readdir(ofdir)
    any(filename -> startswith(filename, iname), dir)
end

function writeMe(iname::String,ivs::DataFrame,A::Array{Float64,2};ofdir=fn.ofDir)
    ofivs= join([iname,nrow(ivs)],"_") * "-init.csv"
    oftrv= join([iname,nrow(ivs)],"_") * "-trav.dat"
    CSV.write(joinpath(ofdir,ofivs),ivs)
    writedlm(joinpath(ofdir,oftrv),A)
end

#locating BTS and OpenFlights input files
const APdir= joinpath("..","data","by-tract","air")::String
#raw BTS data, with separate per-carrier flights
global const ifBTS=joinpath(APdir,"2019 BTS domestic.csv")::String
#raw OpenFlights AP data
global const ifAPs=joinpath(APdir,"Openflights airports.dat")::String
#TODO: consider wgetting from the original https://raw.githubusercontent.com/jpatokal/openflights/master/data/airports.dat

#===FluTE======READ=&=PROCESS===FluTE===TRACTS=============#
# sample usage: Oboe.lsTracts()[4] |> Oboe.rdFluteTract |> Oboe.aggBySte
# a *node* (tract) must have the fields :Pop,:LAT,:LNG

#show the FluTE's tract filenames found in ins.ifDir, default to fn
function lsTracts(ins::NamingSpec = fn)
    filter(s::String -> endswith(s,ins.fltInitSuff),readdir(ins.ifDir))
end

#= load a FluTE census tract info into a DataFrame,
setting the types and colnames 
Default to NW tracts (Oregon + Washington)
=#
function rdFluteTract(ifName::String=lsTracts()[findfirst(s -> startswith(s,"a~NW"),lsTracts())],
            ins::NamingSpec=fn)
    idf=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    rename!(idf, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
    idf.Name = map((s,z,w)-> join([s,z,w],"~"),idf.Ste,idf.Cty,idf.Tra)
    return idf
end

# read the FluTE's all-US-by-tract instance (usa-tracts.dat)
function rdWholeUS(ins::NamingSpec=fn)
    ifName = filter(s -> split(s,"-")[1]=="usa",lsTracts(ins) )[1]
    rdFluteTract(ifName,ins)
end

# given a DataFrame of tracts, remove those that aren't from the given states
function censorFluteTractByFIPS(tracts::DataFrame, fips::Vector{String})
    filter(row -> row.Ste ∈ fips, tracts)
end





#= Select an `id`-generating function
0. if it's a census tract, just write $Ste~$Cty~$Tra
1. if it's a county or a state, write as Integer to match `us-10m.json`
2. if it's a Voronoi cell of an AP, write the AP's IATA_Code
=#
function select_mkid(nms::Array{String})
    if ["Ste","Cty","Tra"] ⊆ nms #node's a census tract
         return r-> join([r.Ste,r.Cty,r.Tra],"~")
    elseif ["Ste","Cty"] ⊆ nms #node's a county, must accomodate `us-10m.json`
         return r-> parse(Int,r.Ste*"000") + parse(Int,r.Cty)
    elseif  ["Ste"] ⊆ nms #node's a state, must accomodate `us-10m.json`
         return r-> parse(Int,r.Ste)
     elseif ["IATA_Code"] ⊆ nms #node is Voronoi cell around an AP
         return r -> r.IATA_Code
     else
        error("Nodes header not recognized. Can't generate an :id without FIPS or IATA_Code.")
     end
 end


 

#= if the first node is sterile (I_i), seed it with 1 infected (Hi, I'm _idempotent_!)
input MUST have [:I_i,:S_i]=#
function infect!(ns::DataFrame)
    if ns.I_i[1] == 0
        ns[1,:S_i]-=1
        ns[1,:I_i]+=1
    end
    return ns
end

#-----INITIAL---VALUES----------------------#

#=
input MUST have [:Pop,:Name,:LAT,:LNG]
input SHOULD have [:IATA_Code], if absent, inserted as "dummy"
out spec: [:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG]
=#
function ns2iv_sterile(ns::DataFrame)
    out = copy(ns)
    rename!(out, :Pop => :N_i) 
    insertcols!(out,1,:id => map(select_mkid(names(out)), eachrow(out)))
    insertcols!(out,4, :S_i => out.N_i)
    insertcols!(out,5, :I_i => zeros(Int64,nrow(out)))
    insertcols!(out,6, :R_i => zeros(Int64,nrow(out)))
    if "IATA_Code" ∉ names(out)
        insertcols!(out,2,:IATA_Code => repeat(["dummy"],nrow(out)))
    end
    select(out,[:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG])
end

#quick delegate for the most common use case
ns2iv(ns::DataFrame) = ns2iv_sterile(ns) |> infect! 

end