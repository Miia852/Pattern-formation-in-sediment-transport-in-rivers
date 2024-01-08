function run_simulation(run_function, solver)
    run_function(solver, false)
end

function makeTable(benchmark::BenchmarkTools.Trial, Nx, Nt, solver)
    str = replace(string(solver), r"\([^)]*\)" => "")
    
    df = DataFrame(Algorithm = str)
    df.Nx .= Nx
    df.Nt .= Nt
    # df.Algorithm .= string(Tsit5())
    df.gctimes .= median(benchmark.gctimes)
    df.time .= median(benchmark.times)
    df.memory .= median(benchmark.memory)
    df.allocs .= median(benchmark.allocs)

    return df
end

function benchmark(run_func, algorithms)
    benchmark_results = Dict()
    combined_df = DataFrame()

    for solver in algorithms
        @info "Running simulation for Nx=$Nx, Nt=$Nt, solver=$solver"
        try
            benchmark = @benchmark run_simulation($run_func, $solver)
            benchmark_results[(Nx, Nt, solver)] = benchmark
            df = makeTable(benchmark, Nx, Nt, solver)
            append!(combined_df, df)
            @info "Time estimate: $(minimum(benchmark).time) seconds, Memory estimate: $(minimum(benchmark).memory) bytes"
        catch e
            @warn "Simulation failed for Nx=$Nx, Nt=$Nt, solver=$solver. Error: $e"
        end
    end

    output_path = joinpath(@__DIR__, "..", "output", "csv", "benchmark.csv")
    CSV.write(output_path, combined_df)
    # println(combined_df)
end

