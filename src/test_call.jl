include("../../dtd_julia/main.jl")
using Cxx, Libdl
using StatsBase, Statistics
using Test

# this is meant to be run from within the build directory
addHeaderDir("./eigen/")
Libdl.dlopen(pwd() * "/build/libdtd_jl.so", Libdl.RTLD_GLOBAL)
cxxinclude(pwd() * "/interface_jl.hpp")

const ngenes=100
const ncells=5
const nsamples=10

x = rand(ngenes, ncells)
g = rand(ngenes)
gp = rand(ngenes)
y = rand(ngenes, nsamples)
c = rand(ncells, nsamples)

# define in the cpp space
ptr=pointer(x)
rows=size(x, 1)
cols=size(x, 2)
matcols=rows*cols
cxx_x = icxx"ptrToVec($ptr, $matcols);"

ptr=pointer(y)
rows=size(y, 1)
cols=size(y, 2)
matcols=rows*cols
cxx_y = icxx"ptrToVec($ptr, $matcols);"

ptr=pointer(c)
rows=size(c, 1)
cols=size(c, 2)
matcols=rows*cols
cxx_c = icxx"ptrToVec($ptr, $matcols);"

ptrgp=pointer(gp)
rows = size(gp,1)
cxx_gp = icxx"ptrToVec($ptrgp, $rows);"

ptrg=pointer(g)
rows = size(g,1)
cxx_g = icxx"ptrToVec($ptrg, $rows);"

# TODO: test that matrix conserves rows & cols...

# XTGX_inv:
res=icxx"xtgxinv($cxx_x, $cxx_g, $ngenes, $ncells);"

res_v = collect(res)
xtgxi_cxx = reshape(collect(res), (ncells, ncells));

@test isapprox(xtgxi_cxx, xtgx_inv(x, g))

# test stats function cov, var
res=icxx"cov($cxx_g, $cxx_gp);" # returns a number
@test isapprox(res, cov(g, gp))
res=icxx"var($cxx_g);" # returns a number
@test isapprox(res, var(g))

# now we need a model to test:
model = icxx"new dtd::models::GoertlerModel(makeGoertlerModel($cxx_x, $cxx_y, $cxx_c, $ngenes, $ncells, $nsamples));"

# loss function:
cxx_fnval = icxx" evalGoertlerModel($model, $cxx_g);"
@test isapprox(cxx_fnval, eval_L(x,y,g,c))

# gradient
cxx_grad= icxx" gradGoertlerModel($model, $cxx_g);"
cxx_grad = collect(cxx_grad)
@test isapprox(collect(cxx_grad), gradient_L(x,y,g,c))

# wrapper fn:
function gradientfun!(gr, g)
    gr .= gradient_L(x,y,g,c)
end
function fn(g)
    return eval_L(x,y,g,c)
end

# step width determination
cxx_sw = icxx" bb_learning_rate($model, $cxx_g);"
jl_sw = learning_rate_BB!(gradientfun!, g);

# TODO: test clampPos / clampNeg, pos_subspace, norm,..

# check that all preconditions are still met:
g_iter = copy(g) # <- working copy of g
@test all(g_iter .> 0)
@test isapprox(collect(cxx_g), g_iter)

# one iteration of fista
lambda=0.01
maxiter=2 # maxiter = 1 does nothing.
learning_rate=0.1
linesearchspeed = 2
cycles = 5
cxx_fval_after_one_iteration = icxx"solveFista($model, $cxx_g, $lambda, $maxiter);"
println("run fista from julia: g = ", g_iter)
# run fista:
println("run julia impl of fista:")
fista_impl!(g_iter, lambda, maxiter, learning_rate, gradientfun!, soft_thresholding, nesterov_factor, positive_subspace, fn, identity, linesearchspeed, cycles, true)
println("done julia impl of fista.")
# compare results:
@test isapprox(g_iter, collect(cxx_g))
@test isapprox(cxx_fval_after_one_iteration, fn(g_iter))
