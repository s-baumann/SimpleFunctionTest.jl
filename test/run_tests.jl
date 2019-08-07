using FunctionTest
using Test

# Run tests
println("Test with Float64s")
@time @test include("basic_tests.jl")
