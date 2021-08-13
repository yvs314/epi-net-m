"""
Submodule responsible for handling input and output, i.e.
reading tract and travel data from disk and outputting travel matrices.
"""
module IO

export fn,
       writeMe,
       checkIfFilesExist,
       rdAPs,
       rdBTS,
       rdTidyWfsByFIPS,
       rdFluteTract,
       rdWholeUS,
       censorFluteTractByFIPS,
       ns2iv

using DataFrames
#CSV: IO FluTE & my, DelimitedFiles: writing the matrices (.dat)
using CSV,DelimitedFiles

#===TYPES AND NAMING CONVENTIONS===#
"""Data type for holding paths and naming conventions for I/O"""
struct NamingSpec #all fields are String, don't say I didn't warn you
    # e.g. sep="-", sSep="_"; $name_$size-$suff
    sep::String #to the right of last $sep is filename suffix, to the left is the instance name
    sSep::String
    # read from ifDir, write to ofDir
    ifDir::String
 #  ifDirAir::String
    ofDir::String
    # what's after instance's name in its Initial Values file name
    fltInitSuff::String
    myInitSuff::String
end

#the naming conventions I am going to use
#as well as input and output directories
const fn=NamingSpec("-","_"
    ,joinpath("..","data","by-tract","flute")
    ,joinpath("..","data","by-tract")
    ,"tracts.dat","init.csv")

#===OUTPUT===#
"""Return `true` if there are files whose names start with `iname` in `odfir`."""
function checkIfFilesExist(iname::String; ofdir=fn.ofDir)
    dir = readdir(ofdir)
    any(filename -> startswith(filename, iname), dir)
end

"""Save the initial values table `ivs` and the travel matrix `A` to `ofdir`."""
function writeMe(iname::String,ivs::DataFrame,A::Array{Float64,2};ofdir=fn.ofDir)
    ofivs= join([iname,nrow(ivs)],"_") * "-init.csv"
    oftrv= join([iname,nrow(ivs)],"_") * "-trav.dat"
    CSV.write(joinpath(ofdir,ofivs),ivs)
    writedlm(joinpath(ofdir,oftrv),A)
end

#===AIRPORTS===#
#locating BTS and OpenFlights input files
const APdir= joinpath("..","data","by-tract","air")::String
#raw BTS data, with separate per-carrier flights
const ifBTS=joinpath(APdir,"2019 BTS domestic.csv")::String
#raw OpenFlights AP data
const ifAPs=joinpath(APdir,"Openflights airports.dat")::String
#TODO: consider wgetting from the original https://raw.githubusercontent.com/jpatokal/openflights/master/data/airports.dat
const apAllColNames = [:ID,:Name,:City,:Country,:IATA_Code,:ICAO_Code,:LAT,:LNG,:Altitude,:Timezone,:Daylight_Savings,:TZ,:Type,:Source]
#reordered in the order of necessity; 3-letter IATA code as ID
const apRetainedColNames=[:IATA_Code,:LAT,:LNG,:Name,:City,:Country]

"""
Read the OpenFlights.org's airports.dat, and retain only the columns
that will be used; put :IATA_Code,:LAT,:LNG first, these are definitely useful.
"""
#CAVEAT: some missing values in :City, some weird values in :Country
#TODO: resolve this caveat
function rdAPs(ifName=ifAPs::String)
    out=CSV.File(ifName,header=false) |> DataFrame
    rename!(out, apAllColNames)
    #retain only the useful columns
    select!(out,apRetainedColNames)
end

#---READ---BTS--FLIGHTS-------------------------#
"""
Read the BTS file, and retain only :1 Passengers, :5 ORIGIN, and :7 DEST;
set the columns to [:ORG,:DST,:PSG] for uniformity.
"""
function rdBTS(ifName=ifBTS::String)
    rawBTS = select(CSV.read(ifName,DataFrame), :1,:5,:7)
    if map( x -> floor(x)==x, rawBTS[:,1]) |> all #if :PSG are Integer
        rawBTS.PASSENGERS=convert(Array{Int64,1},rawBTS.PASSENGERS)
    else #∃ non-integer passengers-per-year entry
        error("CAVEAT: *non-integer* PASSENGERS in Flights input.")
    end
    rename!(rawBTS,[:PSG,:ORG,:DST])
    return select(rawBTS,[:ORG,:DST,:PSG]) #
end

#===COMMUTER FLOW (*-wf-*.dat files)===#
"""List the `-wf-` filenames found in `ins.ifDir`, default to `fn`"""
function lsWf(ins::NamingSpec = fn)
    filter(s::String -> startswith(s,"usa-wf-"),readdir(ins.ifDir))
end

ls_fipsAll = map(s -> split(split(s,".")[1],"-")[3] |> string,lsWf()) 

"""
Read and tidy all "usa-wf-\$fips.dat" for \$fips in fipss; default to NW: Oregon + Washington.
If the tracts argument is provided, commutes referencing tracts that aren't in it
will be discarded.
"""
function rdTidyWfsByFIPS(fipss::Array{String,1}=["41","53"], tracts::DataFrame=DataFrame(),
     ins::NamingSpec=fn)
    #generate FluTE filename by State's FIPS
    wf_by_FIPS(fips::String) = "usa-wf-$fips.dat"
    
    wfs = map(wf_by_FIPS, fipss)
    wfs2 = [CSV.File(ipath, 
        header = false, 
        types =[String,String,String,String,String,String,Int64]) |> 
    DataFrame for ipath in map(f -> joinpath(ins.ifDir,f), wfs)]

    to_keep(r) = r[4] ∈ fipss && r[3]≠r[6]
    #TIDY: retain only the commute between states in `fipss`; this chucks the REFLEXIVE commute
    wfs3 = map(df -> filter(to_keep,eachrow(df) ) |> DataFrame,wfs2)
    wfs4 = (length(wfs3) > 1) ? reduce(vcat,wfs3) : wfs3[1] #add them all together
    #combine the State,County,Tract triples into single columns
    insertcols!(wfs4, 1, 
        :ORG => map((s,z,w)-> join([s,z,w],"~"),wfs4.Column1, wfs4.Column2,wfs4.Column3))
    insertcols!(wfs4, 2, 
        :DST => map((s,z,w)-> join([s,z,w],"~"),wfs4.Column4, wfs4.Column5,wfs4.Column6))
    rename!(wfs4, :Column7 => :CMT) #these are daily commuters between :ORG and :DST
    select!(wfs4,:ORG,:DST,:CMT) #chuck the unnecessary
    #Chuck commute pairs that have tract IDs not present in tracts
    wfs5 = isempty(tracts) ? # check if tracts has been passed
        wfs4 :
        scrubWfs(tracts, wfs4) |> DataFrame
    sort!(wfs5,[:ORG,:DST])

    return wfs5
end

"""
Return `wfs` excluding the commutes for tracts not present in `ns`.
`ns` MUST have [:Name]; `wfs` MUST have [:ORG,:DST,:CMT]
"""
function scrubWfs(ns::DataFrame,wfs::DataFrame)
    allCmtTracts = (wfs.ORG |> unique) ∪ (wfs.DST |> unique)
     #these tracts are not present in FluTE's `us-tracts.dat`
    notPresent = setdiff(allCmtTracts,ns.Name)
    retwfs = filter(r -> r.ORG ∉ notPresent && r.DST ∉ notPresent,eachrow(wfs))
    return retwfs
end

#===FluTE TRACTS (*-tracts.dat files)===#
# sample usage: Oboe.lsTracts()[4] |> Oboe.rdFluteTract |> Oboe.aggBySte
# a *node* (tract) must have the fields :Pop,:LAT,:LNG

"""List the FluTE's tract filenames found in `ins.ifDir`, default to `fn`"""
function lsTracts(ins::NamingSpec = fn)
    filter(s::String -> endswith(s,ins.fltInitSuff),readdir(ins.ifDir))
end

"""
Load a FluTE census tract info into a `DataFrame`, setting the types and colnames.
Default to NW tracts (Oregon + Washington)
"""
function rdFluteTract(ifName::String=lsTracts()[findfirst(s -> startswith(s,"a~NW"),lsTracts())],
            ins::NamingSpec=fn)
    idf=CSV.File(joinpath(ins.ifDir,ifName)
    ,header=false
    ,types=[String,String,String,Int64,Float64,Float64]) |> DataFrame
    rename!(idf, [:Ste,:Cty,:Tra,:Pop,:LAT,:LNG])
    idf.Name = map((s,z,w)-> join([s,z,w],"~"),idf.Ste,idf.Cty,idf.Tra)
    return idf
end

"""Read the FluTE's all-US-by-tract instance (`usa-tracts.dat`)"""
function rdWholeUS(ins::NamingSpec=fn)
    ifName = filter(s -> split(s,"-")[1]=="usa",lsTracts(ins) )[1]
    rdFluteTract(ifName,ins)
end

"""Given a `DataFrame` of tracts, remove those that are not from the given states."""
function censorFluteTractByFIPS(tracts::DataFrame, fips::Vector{String})
    filter(row -> row.Ste ∈ fips, tracts)
end

#===EPIDEMIC SIMULATION AND INITIAL VALUES===#
"""
input MUST have [:Pop,:Name,:LAT,:LNG]
input SHOULD have [:IATA_Code], if absent, inserted as "dummy"
Initialize a `DataFrame` for SIR with the columns 
`[:id,:IATA_Code,:N_i,:S_i,:I_i,:R_i,:Name,:LAT,:LNG]`
"""
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
"""Initialize a `DataFrame` for SIR and infect the first node."""
ns2iv(ns::DataFrame) = ns2iv_sterile(ns) |> infect!

"""
If the first node is sterile (`I_i`), seed it with 1 infected (Hi, I'm _idempotent_!)
input MUST have `[:I_i,:S_i]`
"""
function infect!(ns::DataFrame)
    if ns.I_i[1] == 0
        ns[1,:S_i]-=1
        ns[1,:I_i]+=1
    end
    return ns
end

"""
Select an `id`-generating function
0. if it's a census tract, just write \$Ste~\$Cty~\$Tra
1. if it's a county or a state, write as Integer to match `us-10m.json`
2. if it's a Voronoi cell of an AP, write the AP's IATA_Code
"""
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

end