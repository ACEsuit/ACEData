using PyCall, JSON, ACEbase

datapath = "/Users/ortner/Dropbox/fitsWithLogging/data"

ase = pyimport("ase")
pmt = pyimport("pymatgen.core")
pmt2ase = pyimport("pymatgen.io.ase")
json = pyimport("json")

using ACE1pack, StaticArrays, JuLIP, LinearAlgebra

#Converts a matriz into an SVector
function matrix2svector(M)
    sv = [SVector{size(M)[1]}(M[:,i]) for i in 1:size(M)[2]]
    return sv
end

#Takes one of ("Cu", "Ge", "Li", "Mo", "Ni", "Si") as an element and either "training" or "test"
#Returns an array of JuLIP Atoms with energies and forces
function getData(element, traintest)
    filepath = joinpath(datapath, element, traintest * ".json")
    f = open(filepath, "r")
    data = json.load(f)

    output = Atoms[]
    for d in data
        #reads the structure using pymatgen
        curr_structure = pmt.Structure.from_dict(d["structure"])

        #convert that to an ASE (python) object
        at = pmt2ase.AseAtomsAdaptor.get_atoms(curr_structure)
        at_JuLIP = Atoms(ASEAtoms(at))

        #enrich it with energy and forces
        E = JuLIP.JData(Inf, 0.0, d["outputs"]["energy"])
        atnum = d["num_atoms"]
        ftmp = zeros(3,atnum) 
        for i in 1:atnum
            ftmp[:,i] = d["outputs"]["forces"][i,:]
        end
        F = JuLIP.JData(Inf, 0.0, matrix2svector(ftmp))
        at_JuLIP.data = Dict("energy"=>E, "force"=>F)

        push!(output, at_JuLIP)
    end
    return output
end

# create the folder ZuoEtAl2020
outpath = joinpath(@__DIR__(), "..", "ZuoEtAl2020")


for sym in ["Cu", "Ge", "Li", "Mo", "Ni", "Si"]
    train = getData(sym, "training")
    test = getData(sym, "test")
    JuLIP.write_extxyz(joinpath(outpath, "$(sym)_train.xyz"), identity.(train))
    JuLIP.write_extxyz(joinpath(outpath, "$(sym)_test.xyz"), identity.(test))
end

# tar -czvf ZuoEtAl2020.tar.gz ZuoEtAl2020
# rm -rf ZuoEtAl2020


# # --- test that the xyz have all the data 
# sym = "Cu"
# test = JuLIP.read_extxyz(joinpath(outpath, "$(sym)_test.xyz"))
# for at in test
#     @assert haskey(at.data, "energy")
#     @assert haskey(at.data, "force")
#     @show at.data["energy"]
# end