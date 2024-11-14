using Test
using Hippocampus

@testset "Utils" begin
    markers = [84, 11, 21, 31, 12, 22, 42, 13, 23, 33]
    timestamps = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    trial_markers, trial_timestamps = Hippocampus.reshape_triggers(markers, timestamps)
    @test size(trial_markers) == (3,3)
    @test trial_markers == [11 21 31;12 22 42;13 23 33]
    @test trial_timestamps ≈ [1.0 2.0 3.0;4.0 5.0 6.0;7.0 8.0 9.0]

    # test that we can catch errors
    # wrong poster number on the second trigger of the first trial
    markers = [84, 11, 22, 31, 12, 22, 42, 13, 23, 33]
    @test_throws "Inconsistent cue markers" Hippocampus.reshape_triggers(markers, timestamps)

    # wrong index trigger on the first trigger on the last trial
    markers = [84, 11, 21, 31, 12, 22, 42, 23, 23, 33]
    @test_throws "Inconsistent main markers" Hippocampus.reshape_triggers(markers, timestamps)
end

@testset "Unity" begin
    udata = cd(@__DIR__) do 
        Hippocampus.UnityData("unity_data.csv")
    end
    @test size(udata.triggers) == size(udata.timestamps) == (1,3)
    @test udata.triggers == [16 26 36]
    @test udata.timestamps ≈ [0.04301112 1.04500456 8.226267400000001]
end

@testset "Ripple markers" begin
    markers = UInt16.([84,12,22,32])
    timestamps = [7.0486, 9.1376, 10.1465, 16.457833333333333]
    rpdata = Hippocampus.RippleData(markers,timestamps)
    @test size(rpdata.triggers) == size(rpdata.timestamps)  == (1,3)
    @test rpdata.triggers == permutedims(markers[2:end]) 
    @test rpdata.timestamps ≈ permutedims(timestamps[2:end])
end
