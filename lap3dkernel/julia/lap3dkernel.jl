# play with julia, using laplace 3d kernel eval.
# Barnett 9/18/18.  After Steven Johnson tweaking, 2/20/20.
# manual SIMD version, Luiz M Faria 9/7/20.

#= comment block
=#

# How to run:
# mutlithreading: must start julia with
# JULIA_NUM_THREADS=8 julia

# then, to run:  include("lap3dkernel.jl")

using LinearAlgebra
using Printf
using BenchmarkTools
using Base.Threads
using SIMD # helps when "manual" vectorization is needed
using VectorizationBase: pick_vector_width # use to to determine appropriate size of register vector

# for turbo David Stein versions below...
using LoopVectorization

function lap3dcharge(y,q,x)       # vec over targs
    T = eltype(y)
    nt = size(x,2)
    pot = zeros(T,1,nt)    # note zeros(nt), col vec, fails to add r later
    prefac = T(1/(4*pi))
    for j in eachindex(q)           # loop over srcs
        r2 = sum((x .- y[:,j]).^2,dims=1)      # sq dist
        r = sqrt.(r2)
        pot +=  prefac * q[j] ./ r
    end
    return pot
end

function lap3dcharge_devec(y,q,x)       # unwrap i,j.   C-style coding, SIMD
    T = eltype(y)
    nt = size(x,2)
    ns = size(y,2)
    pot = zeros(T,1,nt)    # note zeros(nt), col vec, fails to add r later
    prefac = 1/(4*pi)
    @inbounds for i = 1:nt      # targs
        for j = 1:ns          # srcs;  using @simd here makes no difference
            #r2ij = sum((x[:,i] - y[:,j]).^2)      # sq dist - terrible!
            r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            rij = sqrt(r2ij)
            pot[i] +=  prefac * q[j] / rij
        end
    end
    return pot
end

function lap3dcharge_devec_par(y,q,x)   # multi-threaded version of above
    T = eltype(y)
    nt = size(x,2)
    pot = zeros(T,nt)    # note zeros(nt), col vec, fails to add r later
    prefac = 1/(4*pi)
    @threads for i in eachindex(pot)      # targs.
        for j in eachindex(q)   # srcs;  using @simd here makes no difference
            @inbounds r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            #rij = sqrt(r2ij)
            @inbounds pot[i] +=  prefac * q[j] / sqrt(r2ij) #rij
        end
    end
    return pot
end

# next three funcs are manual SIMD version by Luiz M Faria...
function lap3dcharge_devec_par_manualSIMD(y,q,x,V=Vec{4,Float64})   
    xt,yt,zt = x[1,:],x[2,:],x[3,:]
    xs,ys,zs = y[1,:],y[2,:],y[3,:]
    lap3dcharge_devec_par_manualSIMD(xs,ys,zs,xt,yt,zt,q,V)
end

function lap3dcharge_devec_par_manualSIMD(xs,ys,zs,xt,yt,zt,q,V)
    T = eltype(xs)
    @assert length(xs) == length(ys) == length(zs)
    @assert length(xt) == length(xt) == length(zt)
    ns,nt = length(xs), length(xt)
    pot   = zeros(T,nt)    # note zeros(nt), col vec, fails to add r later
    @threads for t in 1:length(V):nt
        Xt_vec = vload(V,xt,t)
        Yt_vec = vload(V,yt,t)
        Zt_vec = vload(V,zt,t)
        pot_vec = _inner_loop(Xt_vec,Yt_vec,Zt_vec,t,xs,ys,zs,q)
        vstore(pot_vec,pot,t)
    end
    return pot
end

function _inner_loop(Xt_vec::T,Yt_vec::T,Zt_vec::T,t,xs,ys,zs,q) where {T}
    prefac = 1/(4*pi)
    pot_vec = zero(T)
    for s in eachindex(xs)
        dx = xs[s] - Xt_vec
        dy = ys[s] - Yt_vec
        dz = zs[s] - Zt_vec
        r2 = dx*dx + dy*dy + dz*dz
        pot_vec +=  prefac * q[s] / sqrt(r2) #rij
    end
    return pot_vec
end

function lap3dcharge_tturbo(y,q,x)   # David Stein LoopVec tturbo
    # NB tturbo means:  @turbo threads=true
    T = eltype(y)
    nt = size(x,2)
    pot = zeros(T,nt)
    prefac = 1/(4*pi)
    @tturbo for i in eachindex(pot)      # targs.
        accum = zero(T)
        for j in eachindex(q)
            r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            accum +=  prefac * q[j] / sqrt(r2ij) #rij
        end
        pot[i] = accum
    end
    return pot
end

function lap3dcharge_turbo(y,q,x)   # David Stein LoopVec turbo + threads
    T = eltype(y)
    nt = size(x,2)
    pot = zeros(T,nt)
    prefac = 1/(4*pi)
    @inbounds @threads for i in eachindex(pot)      # targs.
        accum = zero(T)
        @turbo for j in eachindex(q)
            r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            accum +=  prefac * q[j] / sqrt(r2ij) #rij
        end
        pot[i] = accum
    end
    return pot
end




# ================= MAIN ===================================================
# loop over both single and double precision versions...
for T in (Float32,Float64)   # element type (e.g. Float64 or Float32)
    N = pick_vector_width(T) # number of elements of type T that can be stored on the vector register
    V = Vec{N,T} # SIMD.jl type (used just for lap3dcharge_devec_par_manualSIMD)
    
    ns = 10000
    nt = 10000
    x = rand(T,3,nt)
    y = rand(T,3,ns)
    q = randn(T,ns)    # charges
    t = @elapsed lap3dcharge(y,q,x)    # discards return value
    t = @elapsed lap3dcharge(y,q,x)    # discards return value - is already compiled
    check = lap3dcharge(y,q,x) |> sum
    @printf("Result with type %s: \n",T)
    @printf("targ-vec: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)
    # 0.06 Gpair/s, ie 5x slower than devec

    lap3dcharge_devec(y,q,x)   # compile it?
    t = @elapsed lap3dcharge_devec(y,q,x)    # discards return value
    t = @elapsed lap3dcharge_devec(y,q,x)    # discards return value
    check = lap3dcharge_devec(y,q,x) |> sum
    @printf("devec: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)
    # best, single-threaded,  0.33 Gpair/s   (py+numba was 1.3 since multithreaded)

    t = @elapsed lap3dcharge_devec_par(y,q,x)    # discards return value
    t = @elapsed lap3dcharge_devec_par(y,q,x)    # discards return value
    check = lap3dcharge_devec_par(y,q,x) |> sum
    @printf("devec par: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)
    # 1.26 Gpair/s, matches py+numba

    # Luis Faria's tests of Sept 2020:
    if false          # seems to fail after update to 1.6.2, pkg updates
        t = @elapsed lap3dcharge_devec_par_manualSIMD(y,q,x,V)    # discards return value
        t = @elapsed lap3dcharge_devec_par_manualSIMD(y,q,x,V)    # discards return value
        check = lap3dcharge_devec_par_manualSIMD(y,q,x,V) |> sum
        @printf("devec par new: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)
    end

    # David Stein's tests of July 2021:
    t = @elapsed lap3dcharge_turbo(y,q,x)    # discards return value
    t = @elapsed lap3dcharge_turbo(y,q,x)    # discards return value
    check = lap3dcharge_turbo(y,q,x) |> sum
    @printf("LV turbo: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)

    t = @elapsed lap3dcharge_tturbo(y,q,x)    # discards return value
    t = @elapsed lap3dcharge_tturbo(y,q,x)    # discards return value
    check = lap3dcharge_tturbo(y,q,x) |> sum
    @printf("LV tturbo: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)
end

# S. Johnson chat: use benchmark tools.
# @btime lap3dcharge_devec_par($y,$q,$x);
# the $ is a btime thing: means figure out the type of args before benchmarking
# ; suppresses output - but actually doesn't.

# Ryzen2 laptop 5700U (8-core): 20k * 30k test case:
#Result with type Float32: 
#targ-vec: 600000000 src-targ pairs, ans: -104982.625000 
# 	 time 5.57 s 0.108 Gpair/s
#devec: 600000000 src-targ pairs, ans: -104982.671875 
# 	 time 1.46 s 0.411 Gpair/s
#devec par: 600000000 src-targ pairs, ans: -104982.671875 
# 	 time 0.207 s 2.9 Gpair/s
#LV turbo: 600000000 src-targ pairs, ans: -104982.625000 
# 	 time 0.0613 s 9.79 Gpair/s
#LV tturbo: 600000000 src-targ pairs, ans: -104982.703125 
# 	 time 0.0169 s 35.5 Gpair/s
#Result with type Float64: 
#targ-vec: 600000000 src-targ pairs, ans: 1265455.471453 
# 	 time 7.66 s 0.0784 Gpair/s
#devec: 600000000 src-targ pairs, ans: 1265455.471453 
# 	 time 1.91 s 0.314 Gpair/s
#devec par: 600000000 src-targ pairs, ans: 1265455.471453 
# 	 time 0.264 s 2.27 Gpair/s
#LV turbo: 600000000 src-targ pairs, ans: 1265455.471453 
# 	 time 0.103 s 5.82 Gpair/s
#LV tturbo: 600000000 src-targ pairs, ans: 1265455.471453 
# 	 time 0.0652 s 9.21 Gpair/s

# concl: tturbo wins big, esp at single-prec.
