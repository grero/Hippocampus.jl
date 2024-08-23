function reshape_triggers(markers, timestamps)
    # the first marker is a session start; the remaining come in trios
    nn = length(markers)
    if markers[1] == 84
        nn -= 1
        _markers = markers[2:end]
        _timestamps = timestamps[2:end]
    else
        _markers = markers
        _timestamps = timestamps
    end
    rem(nn,3) == 0 || error("Inconsistent number of markers")
    nt = div(nn,3)
    trial_markers = permutedims(reshape(_markers, 3, nt))
    trial_timestamps = permutedims(reshape(_timestamps,3,nt))

    # sanity check; make sure that the last number digit is the same for each trial
    # and that the succession is 1,2,3 or 1,2,4.
    main_marker = floor.(trial_markers/10.0)
    p1 =sum(sum(main_marker .≈ [1.0 2.0 3.0],dims=2).==3)
    p2 =sum(sum(main_marker .≈ [1.0 2.0 4.0],dims=2).==3)

    p1+p2 == nt || error("Inconsistent main markers")

    # check the the minor markers are the same
    minor_marker = trial_markers - 10.0*main_marker
    p3 = sum(sum(minor_marker .== minor_marker[:,1:1],dims=2).==3)
    p3 == nt || error("Inconsistent cue markers")
    trial_markers, trial_timestamps
end

struct Trial
    i::UInt64
end