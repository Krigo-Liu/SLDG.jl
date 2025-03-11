module Test1d

using Test
using Printf
using LinearAlgebra
using FastGaussQuadrature
using Debugger
using SLDG  # Ensure SLDG is available

"""
    run()

Run the SLDG1d test suite.
"""
function run()
    @testset "SLDG1d Tests" begin
        # Define dummy functions for testing
        ax(x) = 1.0
        fun_init(x) = sin(x)
        # For constant advection with speed 1, the exact solution at time t is sin(x - t).
        exact_function(x, t) = sin(x - t)
    
        # Define a list of grid resolutions.
        grid_resolutions = [20, 40, 80, 160, 320]
    
        # Arrays to store the errors.
        L1_errors = Float64[]
        L2_errors = Float64[]
        Linf_errors = Float64[]
    
        for nx in grid_resolutions
            params = SLDG.SimulationParameters(
                nx,            # Number of grid cells.
                1,             # nk: polynomial degree (solution space has nk+1 nodes).
                6,             # N: number of Gaussian quadrature points.
                0.0,           # Left boundary.
                2π,            # Right boundary.
                20.0,          # Final time.
                0.1,           # CFL number.
                11,            # Number of ghost cells.
                fun_init       # Initialization function.
            )
    
            # Run the solver.
            solution, grid, er1, er2, er3 = SLDG.sldg1d(params, ax, fun_init, exact_function)
    
            push!(L1_errors, er1)
            push!(L2_errors, er2)
            push!(Linf_errors, er3)
        end
    
        # Print the LaTeX table header (optional).
        println("\\hline")
    
        # Loop over grid resolutions and print the formatted rows.
        for i in 1:length(grid_resolutions)
            nx = grid_resolutions[i]
            # Format errors in scientific notation.
            L1_str = @sprintf("%8.2E", L1_errors[i])
            L2_str = @sprintf("%8.2E", L2_errors[i])
            Linf_str = @sprintf("%8.2E", Linf_errors[i])
    
            if i == 1
                # For the first row, no convergence rate is available.
                println(@sprintf("%5d & %10s & & %10s & & %10s & \\\\ \\hline", nx, L1_str, L2_str, Linf_str))
            else
                # Compute the convergence rates using a logarithmic ratio.
                rate_L1 = log(L1_errors[i-1] / L1_errors[i]) / log(2)
                rate_L2 = log(L2_errors[i-1] / L2_errors[i]) / log(2)
                rate_Linf = log(Linf_errors[i-1] / Linf_errors[i]) / log(2)
                println(@sprintf("%5d & %10s & %8.2f & %10s & %8.2f & %10s & %8.2f \\\\ \\hline", 
                    nx, L1_str, rate_L1, L2_str, rate_L2, Linf_str, rate_Linf))
            end
        end
    end
end

# When this file is executed directly, run the test.
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running SLDG 1d test")
    run()
end

end  # module Test1d
