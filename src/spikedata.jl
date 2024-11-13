struct Spiketrain
    timestamps::Vector{Float64}
    reference::String
    components::Int64
end

function Spiketrain(fname::String)
    q = MAT.matread(fname)
    Spiketrain(q["timestamps"][:], q["reference"], q["components"])
end

function Spiketrain()
    if isfile("spiketrain.mat")
        return Spiketrain("spiketrain.mat")
    elseif isfile("spiketrain.csv")
        timestamps = parse.(Float64, readlines(open("spiketrain.csv")))
        return Spiketrain(timestamps, "", 0)
    else
        return Spiketrain(Float64[], "",0)
    end
end

function Base.show(io::IO, x::Spiketrain)
    nspikes = length(x.timestamps)
    print(io, "Spiketrain with $(nspikes) spikes")
end