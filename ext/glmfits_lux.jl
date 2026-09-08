module GLMFitLux
using RecurrentNetworkModels
using RecurrentNetworkModels: Lux, Reactant
using RecurrentNetworkModels.Lux: Training
using RecurrentNetworkModels.Optimisers: Adam
using RecurrentNetworkModels:MLDataDevices

using Hippocampus

function Hippocampus.lossfunc_lux(model, ps, st, (x,y,L))
    η,st_new = model(x, ps, st)
    ll = sum(RecurrentNetworkModels.poisson_loss.(η,y))
    @show typeof(η) typeof(L)
    ll2 = 0.01f0*η*L*η'
    ll = ll + ll2
    return ll, st_new, (;y_pred=η)
end 

function Hippocampus.lossfunc_lux(η, y)
    ll = sum(RecurrentNetworkModels.poisson_loss.(η,y))
    return ll
end

"""
GLM fit using Lux
"""
function Hippocampus.fit_glm_lux(X::Matrix{Float32}, y::Vector{<:Integer},L::Matrix{Float32};nepochs=100,learning_rate=Float32(1e-3))

    dev = MLDataDevices.reactant_device()
    cdev = MLDataDevices.cpu_device()

    d,n = size(X)
    model = Lux.Dense(d,1)
    rng = Random.default_rng()
    _ps,_st = Lux.setup(rng, model)
    ps,st = dev((_ps, _st))
    train_state = Training.TrainState(model, ps, st, Adam(learning_rate))

    #mini batch
 

    # validation set
    nval = round(Int64, 0.2*n) 
    validx = shuffle(1:n)[1:nval]
    trainidx = setdiff(1:n, validx)
    xv = X[:,validx]
    yv = reshape(y[validx], 1,nval)
    (xe,ye,Le) = dev.((xv,yv,L[validx,validx]))
    model_compiled = Reactant.@compile model(xe, ps, Lux.testmode(st))

    sort!(validx)
    nt = round(Int64, 0.5*(n-nval))
    ll = lossfunc_lux(model_compiled, ps, st, (xe,ye,Le))
    @show ll[1]
    prog = Progress(nepochs, "Training...")
    #lossfunc = RecurrentNetworkModels.LPoissonLoss()
    for i in 1:nepochs
        tidx = shuffle(trainidx)[1:nt]
        xt = X[:,tidx]
        yt = reshape(y[tidx], 1, length(tidx))
        (xtd,ytd,Ld) = dev.((xt, yt,L[tidx,tidx]))
         (_, loss, _, train_state) = Training.single_train_step!(
                Training.AutoEnzyme(), lossfunc, (xtd, ytd,Ld), train_state
            )
        st_test= Lux.testmode(train_state.states) #
        η_,_ = model_compiled(xe, train_state.parameters, st_test)
        η = cdev(η_)
        lval = lossfunc(η,yv)
        next!(prog;showvalues=[(:lval, lval)])
    end
    model, cdev((train_state.parameters, train_state.states))
end

end