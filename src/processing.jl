using DataFrames
using CSV

hasempty(row) = "Empty" in row || "Emptyrow" in row || "EmptyRow" in row
matchfile(;kwargs...) = row -> all(getproperty(row, k) == v for (k, v) in kwargs)
const index = filter(r -> !hasempty(r), CSV.read(DATADIR * "file_description.csv", DataFrame))

function search(strain, backbone, plasmid, iptg)
    f = matchfile(strain=strain, backbone=backbone, plasmid=plasmid, iptg=iptg)
    first(filter(f, index))
end

function search(strain, backbone, plasmid)
    f = matchfile(strain=strain, backbone=backbone, plasmid=plasmid)
    sort(filter(f, index), [:iptg])
end

function search(;kwargs...)
    f = matchfile(;kwargs...)
    res = filter(f, index)
    first(size(res)) == 1 ? first(res) : res
end
