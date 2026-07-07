using Pkg
Pkg.activate(".")
using LLRParsing
using DelimitedFiles

function copy_info_files(metadata_file)

    metadata = readdlm(metadata_file, ',', String, skipstart = 1)

    for relpath in metadata[:,1]

        infofile = joinpath(relpath,"base","info.csv")

        path = "tmp/su4/info"
        name = basename(relpath)
        ispath(joinpath(path)) || mkpath(joinpath(path))
        cp(infofile,joinpath(path,"$(name)_info.csv"))
    end
end

metadata_file = "metadata/runs_su4.csv"
copy_info_files(metadata_file)
