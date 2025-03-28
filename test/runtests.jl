using Test
using Hippocampus
using StatsBase
using StableRNGs

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

    # test disk
    f = Hippocampus.disk(1)
    @test f ==  [0.0 0.2 0.0; 0.2 0.2 0.2; 0.0 0.2 0.0]
end

@testset "Paths" begin
    results = Hippocampus.find_path(6, 1,2;path_length=22, max_npaths=1)
    @test length(results) == 1
    @test results[1][1] == 1
    @test results[1][end] == 2
    @test length(results[1]) == 23
    edges = [(i,j) for (i,j) in zip(results[1][1:end-1], results[1][2:end])]
    edge_count = countmap(edges)
    @test unique(values(edge_count)) == [1]

    # test default parameters
    results = Hippocampus.find_path(6,1,2)
    @test length(results) == 1
    @test length(results[1]) == 29 
end

@testset "Plot functions" begin
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
    @test Hippocampus.numtrials(udata) == 1

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
    @test length(bins[1][3]) == 51
    @test bins[1][3][[1,51]] ≈ [0.0,5.0]

    @test bins[2][1][[1,51]] ≈ [-7.5,-2.5]
    @test bins[2][2] ≈ [7.4,7.6]
    @test length(bins[1][3]) == 51
    @test bins[2][3][[1,51]] ≈ [0.0,5.0]

    @test bins[3][1] ≈ [-7.6,-7.4]
    @test bins[3][2][[1,51]] ≈ [2.5,7.5]
    @test length(bins[3][3]) == 51
    @test bins[3][3][[1,51]] ≈ [0.0, 5.0]

    @test bins[4][1] ≈ [-2.6,-2.4]
    @test bins[4][2][[1,51]] ≈ [2.5,7.5]
    @test length(bins[3][3]) == 51
    @test bins[3][3][[1,51]] ≈ [0.0, 5.0]

    #poster assignment
    _poster_pos = Dict(:camel => (-5.0, 1.5, -7.55), :cat => (-7.55, 1.5, 5.0), :rabbit => (-5.0, 1.5, -2.45), :donkey => (5.0, 1.5, 7.55), :croc => (7.55, 1.5, -5.0), :pig => (5.0, 1.5, 2.45))
    mm = Hippocampus.MazeModel()
    posters = Hippocampus.Posters(mm,_poster_pos)
    _pos = mean(posters.sprite[1].points)
    @test _pos[1] ≈ -5.0
    @test _pos[2] ≈ -7.55
    @test _pos[3] ≈ 2.4999998

    #visualizing the maze
    fig = Hippocampus.visualize([(mm,posters)];show_normals=true)
    # make sure we are getting a label on top and a single scene below that
    @test length(fig.content) == 2
    @test typeof(fig.content[1]) <: Hippocampus.Makie.Label
    @test typeof(fig.content[2]) <: Hippocampus.Makie.LScene
    @test length(fig.content[2].scene.plots) == 51 
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

@testset "Spatial" begin
    rpdata = cd(@__DIR__) do
        Hippocampus.RippleData()
    end
    udata = cd(@__DIR__) do
        Hippocampus.UnityData()
    end
    # spatial occupancy
    xbins = range(-12.5, stop=12.5, length=40);
    ybins = xbins
    spoc = Hippocampus.SpatialOccupancy(udata, xbins, ybins)
    @test size(spoc.weight) == (39,39)
    ee = Hippocampus.compute_entropy(spoc)
    @test ee ≈ 8.336202756717388

    # test consistency
    @test Hippocampus.numtrials(rpdata) == Hippocampus.numtrials(udata) == 63
    @test rpdata.triggers == udata.triggers

    # create a spike train
    rng = StableRNG(1234)
    spikes = Hippocampus.model_place_field(udata, rpdata;rng=rng)
    @test length(spikes) == 216

    # create a spatial representation
    sptrain = Hippocampus.Spiketrain(1000.0*spikes, "",0)
    spr = Hippocampus.SpatialRepresentation(sptrain, rpdata, udata)
    # test some random position
    @test spr.position[5][1] ≈ [2.0808, -0.5886]

    # create a spatial map
    spm = Hippocampus.SpatialMap(spr, xbins,ybins, spoc)
    ee = Hippocampus.compute_entropy(spm)
    @test ee ≈ 6.361282290710195
    sic = Hippocampus.compute_sic(spm)
    # TODO: Is this value actually accurate?
    @test sic ≈ 2.103883351775456
end

@testset "Raytrace" begin
    gx = 950.4f0
    gy = 1069.1f0
    # center and normalize
    gx -= 1919.0/2
    gx /= 1919.0
    gy -= 1079/2
    gy /= 1079.0

    # position in the maze
    px = 0.0
    py = -10.0

    gmx,gmy,gmz = Hippocampus.raytrace(gx,gy,[px,py],0.0,60.0, 0.3)

    @test gmx ≈ -0.060816102664816116 
    @test gmy ≈ -3.760296369084596
    @test gmz ≈ 4.926764210409394 
end

@testset "Gaze histogram" begin
    # create a random gaze path on the inner surfaces of the maze
    # and check that we can recover that path when computing histograms
    points = fill(0.0, 3, 10)
    points[:,1] .= [2.5, 7.51, 1.5] # point on the red pillar near the donkey poster
    points[:,2] .= [7.51, -2.5, 1.5] # point on the green pillar near the rabiit poster
    points[:,3] .= [-2.5, -7.51, 1.5] # point on the blue pillar near the cat poster
    points[:,4] .= [-7.51, 2.5, 1.5] # point  on the yellow pillar near the camel poster
    points[:,5] .= [-2.5, -2.49, 1.5] # point on the blue pillar near the pig poster
    points[:,6] .= [2.5, 2.49, 1.5] # point on the red pillar near the crocodile poster
    points[:,7] .= [0.0, 0.0, 4.93] # point in the middle of the ceiling
    points[:,8] .= [-2.72, 12.63, 1.5] # point on one of the peripheral walls
    points[:,9] .= [-2.72, -12.37, 1.5] # point on another peripheral wall
    points[:,10] .=[0.0, 0.0, 0.0] # point on the middle of the floor

    # create the maze
    mm = Hippocampus.MazeModel(;xmin=-12.72, xmax=12.28,ymin=-12.37, ymax=12.63,Δ=0.05,n_vertical_wall_bins=8, n_vertical_pillar_bins=5,
                                n_horizontal_pillar_bins=8)
    nbins = Hippocampus.num_bins(mm)
    @test nbins == 5186
    bins = Hippocampus.get_bins(mm)
    counts,_idx = Hippocampus.compute_histogram([points], bins) 
    # unpack the trial again
    idx = _idx[1]
    # test that the labels are correct
    @test idx[1][4] == :pillars
    @test idx[2][4] == :pillars
    @test idx[3][4] == :pillars
    @test idx[4][4] == :pillars
    @test idx[4][4] == :pillars
    @test idx[6][4] == :pillars
    @test idx[7][4] == :ceiling
    @test idx[8][4] == :walls
    @test idx[9][4] == :walls
    @test idx[10][4] == :floor

end

@testset "Path on cube" begin
    r1 = Rect3f(0.0, 0.0, 0.0,1.0, 1.0, 1.0)
    d = Hippocampus.distance(Point3f(0.5, 0.0, 0.5), Point3f(0.5, 1.0, 0.5), r1)
    @test d ≈ 2.0f0
end