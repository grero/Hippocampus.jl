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
    # test tuple parsing
    pos = parse(NTuple{3,Float64}, "(1.0, 1.2, 3.0)")
    @test pos == (1.0, 1.2, 3.0)
    # test position parsing
    j,pos = Hippocampus.parse_poster_location("P1(-5,1.5,-7.55)")
    @test j == 1
    @test pos == (-5.0,1.5,-7.55)

    # test poster index
    pidx = Hippocampus.get_poster_index(" DTL Donk")
    @test pidx == :donkey 

    udata = cd(@__DIR__) do 
        Hippocampus.UnityData("unity_data.csv")
    end
    @test size(udata.triggers) == size(udata.timestamps) == (1,3)
    @test udata.triggers == [16 26 36]
    @test udata.timestamps ≈ [0.04301112 1.04500456 8.226267400000001]

    tt,pp1,pp2, hh = Hippocampus.get_trial(udata, 1)
    ll = length(tt)
    @test ll == 207
    @test tt[[1,ll]] ≈ [0.04301112, 8.226267400000001]
    @test pp2[[1,ll]] ≈ [-10.0, -1.5189]
    @test pp1[[1,ll]] ≈ [0.0, -6.0951]
    @test hh[[1,ll]] ≈ [0.0, 171.2082]

    # test soft_range
    rr = Hippocampus.soft_range(0.0, 3.1, 1.1)
    @test step(rr) ≈ 1.0333333333333334

    # maze building tests
    bins,normals = Hippocampus.create_mesh((-7.5, 2.5), (-2.5, 7.5), 0.1)
    @test normals ≈ [[0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    @test bins[1][1][[1,51]]  ≈ [-7.5, -2.5]
    @test bins[1][2] ≈ [2.4,2.6]
    @test bins[1][3][[1,10]] ≈ [0.0,5.0]

    @test bins[2][1][[1,51]] ≈ [-7.5,-2.5]
    @test bins[2][2] ≈ [7.4,7.6]
    @test bins[2][3][[1,10]] ≈ [0.0,5.0]

    @test bins[3][1] ≈ [-7.6,-7.4]
    @test bins[3][2][[1,51]] ≈ [2.5,7.5]
    @test bins[3][3][[1,10]] ≈ [0.0, 5.0]

    @test bins[4][1] ≈ [-2.6,-2.4]
    @test bins[4][2][[1,51]] ≈ [2.5,7.5]
    @test bins[3][3][[1,10]] ≈ [0.0, 5.0]
end

@testset "Ripple markers" begin
    markers = UInt16.([84,12,22,32])
    timestamps = [7.0486, 9.1376, 10.1465, 16.457833333333333]
    rpdata = Hippocampus.RippleData(markers,timestamps)
    @test size(rpdata.triggers) == size(rpdata.timestamps)  == (1,3)
    @test rpdata.triggers == permutedims(markers[2:end]) 
    @test rpdata.timestamps ≈ permutedims(timestamps[2:end])
end

@testset "Eyelink" begin
    ee0 = Hippocampus.zerounless(Hippocampus.Eyelink.Event,message="Start Trial 11",sttime=zero(UInt32))
    @test ee0.message == "Start Trial 11"
    @test ee0.sttime == zero(UInt32)
    ee1 = Hippocampus.zerounless(Hippocampus.Eyelink.Event,message="Cue Offset 21",sttime=UInt32(1000))
    ee2 = Hippocampus.zerounless(Hippocampus.Eyelink.Event,message="End Trial 31",sttime=UInt32(50000))
    messages = [ee0,ee1,ee2]
    trial_markers, trial_timestamps = Hippocampus.get_markers(messages)
end
