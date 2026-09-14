# Fixed mean degree and spells per actor: 4× vertices should remain below
# quadratic (16×) runtime/memory growth. Allow headroom for heap depth and
# random reachability variation; keep these checks separate from absolute time.
using Test
include("benchmarks.jl")
results = run(SUITE; verbose=false, samples=30, evals=1, seconds=1)
small = median(results["earliest_arrival"]["sweep_n500"])
large = median(results["earliest_arrival"]["sweep_n2000"])
time_ratio = BenchmarkTools.time(large) / BenchmarkTools.time(small)
memory_ratio = BenchmarkTools.memory(large) / BenchmarkTools.memory(small)
println("SCALING\tearliest_arrival\t4x_vertices\ttime_ratio=", time_ratio,
        "\tmemory_ratio=", memory_ratio)
@testset "Temporal path scaling" begin
    @test time_ratio < 12
    @test memory_ratio < 8
end
