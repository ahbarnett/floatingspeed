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
    
ns = 10000
nt = 10000
x = rand(3,nt)
y = rand(3,ns)
q = randn(ns)    # charges
t = @elapsed lap3dcharge(y,q,x)    # discards return value 
t = @elapsed lap3dcharge(y,q,x)    # discards return value - is already compiled
#pot, t = @timed lap3dcharge(y,q,x)      # gets the time in secs
@printf("targ-vec: %d src-targ pairs in %.3g s: %.3g Gpair/s\n",ns*nt,t,ns*nt/t/1e9)
# 0.06 Gpair/s, ie 5x slower than devec

lap3dcharge_devec(y,q,x)   # compile it?
t = @elapsed lap3dcharge_devec(y,q,x)    # discards return value
t = @elapsed lap3dcharge_devec(y,q,x)    # discards return value
@printf("devec: %d src-targ pairs in %.3g s: %.3g Gpair/s\n",ns*nt,t,ns*nt/t/1e9)
# best, single-threaded,  0.33 Gpair/s   (py+numba was 1.3 since multithreaded)

t = @elapsed lap3dcharge_devec_par(y,q,x)    # discards return value
t = @elapsed lap3dcharge_devec_par(y,q,x)    # discards return value
@printf("devec par: %d src-targ pairs in %.3g s: %.3g Gpair/s\n",ns*nt,t,ns*nt/t/1e9)
# 1.26 Gpair/s, matches py+numba

# S. Johnson chat: use benchmark tools.
# @btime lap3dcharge_devec_par($y,$q,$x);
# the $ is a btime thing: means figure out the type of args before benchmarking
# ; suppresses output - but actually doesn't.
