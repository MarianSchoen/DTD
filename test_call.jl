include("../dtd_julia/main.jl")
using Cxx, Libdl
using StatsBase, Statistics
using Test

# this is meant to be run from within the build directory
addHeaderDir("./src/")
addHeaderDir("/usr/include/eigen3/")
Libdl.dlopen(pwd() * "/build/libdtd_jl.so", Libdl.RTLD_GLOBAL)
cxxinclude(pwd() * "/src/interface_jl.hpp")

const ngenes=10
const ncells=5

x = rand(ngenes, ncells)
g = rand(ngenes)
print(x)

# TODO: test that matrix conserves rows & cols...


# function jlmat(a :: Array{Float64,2})
a = x
ptr=pointer(a)
rows=size(a, 1)
cols=size(a, 2)
matcols=rows*cols
# cxx_v = icxx"JuliaMatrix($ptr, $rows, $cols);"
cxx_x = icxx"ptrToVec($ptr, $matcols);"
println("from julia: x=", collect(cxx_x));
println("from julia: x=", collect(cxx_x));
ptrg=pointer(g)
rows = size(g,1)
cxx_g = icxx"ptrToVec($ptrg, $rows);"
println("from julia: g=", collect(cxx_g));

res=icxx"xtgxinv($cxx_x, $cxx_g, $ngenes, $ncells);"
res_v = collect(res)
xtgxi_cxx = reshape(collect(res), (ncells, ncells));
println("from julia: res: ", xtgxi_cxx);

