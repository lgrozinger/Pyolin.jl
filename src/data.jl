using CSV
using DataFrames
using ZipFile
using FCSFiles
using FileIO

DATADIR = "/home/lewis/data/pyolin_dataset/"
DF = CSV.read("file_description.csv", DataFrame)

function filenamefor(;kwargs...)
    DF[reduce(.&, [getproperty(DF, key) .== value for (key, value) in kwargs]), :filename]
end

plt = plot()
for iptg in [0, 10, 100, 500, 1000]
    fn = filenamefor(strain="KT2440", plasmid="1717", backbone="pSeva221", iptg=iptg)[1]
    run = load(DATADIR * fn)
    h = histogram(run["B1-H"]; normalize=true, bins=0:50)
    plot!(plt, h.series_list[2].plotattributes[:y][1:50]; label="$(iptg)")
end



