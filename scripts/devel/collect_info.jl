using Pkg
Pkg.activate(".")
using LLRParsing
using DelimitedFiles

function copy_info_files(metadata_file,path,file)

    metadata = readdlm(metadata_file, ',', String, skipstart = 1)
    io_combined =  open(file,"w")

    for relpath in metadata[:,1]
        infofile = joinpath(relpath,"base","info.csv")
        name = basename(relpath)
        ispath(joinpath(path)) || mkpath(joinpath(path))
        # copy to new path
        cp(infofile,joinpath(path,"$(name)_info.csv"))
        # create a combined file with info on all runs
        for l in eachline(infofile)
            startswith(l,'#') || println(io_combined,l)
        end
    end
    close(io_combined)
    # now remove duplicates to only have a single header 
    lines = unique(readlines(file))
    io_combined =  open(file,"w")
    for l in lines
        println(io_combined,l)
    end
    close(io_combined)
    # now sort the data by volume and replicas
    data, header = readdlm(file,',',header=true)
    data = sortslices(data; dims = 1, by = row -> (row[2],row[3],row[4]))
    io_combined =  open(file,"w")
    writedlm(io_combined,header,',')
    writedlm(io_combined,data,',')
    close(io_combined)
end

path = "tmp/su4/info"
file = "tmp/su4/info_combined.csv"
metadata_file = "metadata/runs_su4.csv"
copy_info_files(metadata_file,path,file)
