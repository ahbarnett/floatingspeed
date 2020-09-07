# play with julia, using laplace 3d kernel eval.
# Barnett 9/18/18.  After Steven Johnson tweaking, 2/20/20.

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
using SIMD
# this is useless, apparently:
#using Devectorize

function eye(n)
    return Matrix(1.0I,n,n)
end

function lap3dcharge(y,q,x)       # vec over targs
    nt = size(x,2)
    pot = zeros(1,nt)    # note zeros(nt), col vec, fails to add r later
    prefac = 1/(4*pi)
    for j in eachindex(q)           # loop over srcs
        r2 = sum((x .- y[:,j]).^2,dims=1)      # sq dist
        r = sqrt.(r2)
        pot +=  prefac * q[j] ./ r
    end
    return pot
end

function lap3dcharge_devec(y,q,x)       # unwrap i,j.   C-style coding, SIMD
    nt = size(x,2)
    ns = size(y,2)
    pot = zeros(1,nt)    # note zeros(nt), col vec, fails to add r later
    prefac = 1/(4*pi)
@inbounds   for i = 1:nt      # targs
        @simd for j = 1:ns          # srcs
            #r2ij = sum((x[:,i] - y[:,j]).^2)      # sq dist - terrible!
            r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            rij = sqrt(r2ij)
            pot[i] +=  prefac * q[j] / rij
        end
    end
    return pot
end

function lap3dcharge_devec_par(y,q,x)   # multi-threaded version of above
    nt = size(x,2)
    pot = zeros(nt)    # note zeros(nt), col vec, fails to add r later
    prefac = 1/(4*pi)
    @threads for i in eachindex(pot)      # targs.
        @simd for j in eachindex(q)          # srcs
            @inbounds r2ij = (x[1,i]-y[1,j])^2+ (x[2,i]-y[2,j])^2+ (x[3,i]-y[3,j])^2
            #rij = sqrt(r2ij)
            @inbounds pot[i] +=  prefac * q[j] / sqrt(r2ij) #rij
        end
    end
    return pot
end

function lap3dcharge_devec_par_new(y,q,x,V=Vec{4,Float64})   # multi-threaded version of above
    xt,yt,zt = x[1,:],x[2,:],x[3,:]
    xs,ys,zs = y[1,:],y[2,:],y[3,:]
    lap3dcharge_devec_par_new(xs,ys,zs,xt,yt,zt,q,V)
end

function lap3dcharge_devec_par_new(xs,ys,zs,xt,yt,zt,q,V)
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
    
ns = 10000
nt = 10000
x = rand(3,nt)
y = rand(3,ns)
q = randn(ns)    # charges
t = @elapsed lap3dcharge(y,q,x)    # discards return value 
t = @elapsed lap3dcharge(y,q,x)    # discards return value - is already compiled
#pot, t = @timed lap3dcharge(y,q,x)      # gets the time in secs
check = lap3dcharge(y,q,x) |> sum
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

t = @elapsed lap3dcharge_devec_par_new(y,q,x)    # discards return value
t = @elapsed lap3dcharge_devec_par_new(y,q,x)    # discards return value
check = lap3dcharge_devec_par_new(y,q,x) |> sum
@printf("devec par new: %d src-targ pairs, ans: %f \n \t time %.3g s %.3g Gpair/s\n",ns*nt,check,t,ns*nt/t/1e9)


# S. Johnson chat: use benchmark tools.
# @btime lap3dcharge_devec_par($y,$q,$x);
# the $ is a btime thing: means figure out the type of args before benchmarking
# ; suppresses output - but actually doesn't.
