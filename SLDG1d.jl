module SLDG1d

using LinearAlgebra
using FastGaussQuadrature
using Debugger


# -----------------------------
# Define custom types
# -----------------------------
mutable struct Vertex
    coor::Float64
    id::Int
end

mutable struct EulerianElement
    umodal::Vector{Float64}  # size: nk+1
    xgl::Vector{Float64}       # Gauss-Lobatto nodes (size: nk+1)
end

mutable struct Segment
    porigin::Vertex
    pend::Vertex
    id::Int
end

mutable struct UpstreamElement
    umodal::Vector{Float64}   # size: nk+1 (to be computed)
    xgl_star::Vector{Float64} # Gauss-Lobatto nodes at upstream locations (size: nk+1)
    point_origin::Vertex
    point_end::Vertex
    point_inter::Vector{Vertex}  # length: mx+2
    segment::Vector{Segment}     # segments in the element
    nsub::Int
end

# -----------------------------
# Simulation configuration (inputs)
# -----------------------------
struct SimulationParameters
    nx::Int64           # Number of grid points
    nk::Int64           # Polynomial degree (solution space has nk+1 nodes)
    N::Int64            # Number of Gaussian quadrature points
    xleft::Float64      # Left boundary
    xright::Float64     # Right boundary
    time_final::Float64 # Final time
    cfl::Float64        # CFL number
    numberOfGhostCells::Int64  # Number of ghost cells
    InitializationFunction::Function  # Function to initialize the solution
end

# -----------------------------
# Simulation state (internal evolving state)
# -----------------------------
mutable struct SimulationState
    nx::Int64                       
    nk::Int64                       
    numberOfGhostCells::Int64       
    N::Int64                        
    dx::Float64                     
    pi::Float64                     
    eps::Float64                    
    coefficientsForCorrection::Vector{Float64} 
    xleft::Float64                  
    xright::Float64                 
    dt::Float64                     
    time::Float64                   
    time_final::Float64             
    cfl::Float64                    
    x::Vector{Float64}              
    xGrid::Vector{Float64}          
    vertices::Vector{Vertex}                
    vertex_star::Vector{Vertex}             
    eulerianElements::Vector{EulerianElement}   
    upstreamElements::Vector{UpstreamElement}     
end

# -----------------------------
# Initialize fundamental constants and quadrature data.
# -----------------------------
function parameters!(sim::SimulationState)
    sim.pi = 4 * atan(1.0)
    sim.eps = eps(Float64) * 10
    # Set coefficients for correction (assuming nk <= 4)
    if sim.nk == 1
        sim.coefficientsForCorrection .= [1.0, 12.0]
    elseif sim.nk == 2
        sim.coefficientsForCorrection .= [1.0, 12.0, 180.0]
    elseif sim.nk == 3
        sim.coefficientsForCorrection .= [1.0, 12.0, 180.0, 2800.0]
    elseif sim.nk == 4
        sim.coefficientsForCorrection .= [1.0, 12.0, 180.0, 2800.0, 44100.0]
    else
        error("Unsupported nk value")
    end
end

# -----------------------------
# Set up the simulation state (grid, domain, and allocate arrays).
# -----------------------------
function setup!(sim::SimulationState, parameters::SimulationParameters)
    sim.nx = parameters.nx
    sim.nk = parameters.nk
    sim.numberOfGhostCells = parameters.numberOfGhostCells
    sim.xleft = parameters.xleft
    sim.xright = parameters.xright
    sim.cfl = parameters.cfl
    sim.time_final = parameters.time_final
    sim.N = parameters.N

    # Compute the spatial step size.
    sim.dx = (sim.xright - sim.xleft) / sim.nx

    # Allocate basic grid arrays.
    sim.x = [sim.xleft + (i - 0.5)*sim.dx for i in (1-sim.numberOfGhostCells):(sim.nx+sim.numberOfGhostCells)]
    sim.xGrid = [sim.xleft + (i - 1)*sim.dx for i in (1-sim.numberOfGhostCells):(sim.nx+1+sim.numberOfGhostCells)]

    # Allocate vertices (cell boundaries) as Vertex objects.
    sim.vertices = [Vertex(sim.xleft + (i - 1)*sim.dx, i) for i in 1:(sim.nx+1)]
    # Initialize vertex_star as a copy of vertices.
    sim.vertex_star = [Vertex(v.coor, v.id) for v in sim.vertices]

    # Allocate Eulerian elements: one per cell.
    sim.eulerianElements = [EulerianElement(zeros(sim.nk+1), zeros(sim.nk+1)) for _ in 1:sim.nx+sim.numberOfGhostCells*2]
    for i in 1+sim.numberOfGhostCells:sim.nx+sim.numberOfGhostCells
        for j in 1:(sim.nk+1)
            sim.eulerianElements[i].xgl[j] = gausslobatto(sim.nk+1)[1][j]/2*sim.dx + sim.x[i]
        end
    end

    # Allocate Upstream elements: one per cell.
    sim.upstreamElements = [UpstreamElement(zeros(sim.nk+1), zeros(sim.nk+1), Vertex(0.0, 0), Vertex(0.0, 0), Vertex[], Segment[], 0) for _ in 1:sim.nx+sim.numberOfGhostCells*2]

    # Initialize time variables.
    sim.dt = 0.0
    sim.time = 0.0
end

# -----------------------------
# Initialize the simulation (set initial conditions and element values)
# -----------------------------
function initialize!(sim::SimulationState, parameters::SimulationParameters, InitializationFunction::Function)
    # Initialize Eulerian elements using the InitializationFunction.

    for i in 1+sim.numberOfGhostCells:sim.nx+sim.numberOfGhostCells
        for kk in 0:sim.nk
            utemp = 0.0
            for ig in 1:(sim.nk+1)
                # utemp += InitializationFunction(sim.x[i + sim.numberOfGhostCells]+ (gausslegendre(sim.N))[1][ig]/2 * sim.dx) * fle(kk, (gausslegendre(sim.N))[1][ig]/2) * (gausslegendre(sim.N))[2][ig] / 2
                utemp += InitializationFunction(sim.x[i ]+ (gausslegendre(sim.nk+1))[1][ig]/2 * sim.dx) * fle(kk, (gausslegendre(sim.nk+1))[1][ig]/2) * (gausslegendre(sim.nk+1))[2][ig] / 2
            end
           sim.eulerianElements[i].umodal[kk+1] = utemp * sim.coefficientsForCorrection[kk+1]
        end
    end
end

# -----------------------------
# Update time step dt based on the CFL condition and update current time.
# -----------------------------
function update_dt(sim::SimulationState)
    sim.dt = sim.cfl * sim.dx
    if sim.time + sim.dt > sim.time_final
        sim.dt = sim.time_final - sim.time
    end
    sim.time += sim.dt
end

# -----------------------------
# Perform a Runge–Kutta step.
# -----------------------------
function runge_kutta(vx::Float64, dt::Float64, ax::Function)
    vx1 = vx - 0.5 * dt * ax(vx)
    vx2 = vx - 0.5 * dt * ax(vx1)
    vx3 = vx - dt * ax(vx2)
    vx_star = (-vx + vx1 + 2*vx2 + vx3) / 3 - (ax(vx3) * dt) / 6
    return vx_star
end

# -----------------------------
# Apply boundary conditions to the Eulerian elements.
# -----------------------------
function boundary!(sim::SimulationState)
    temp_eulerianElements = [EulerianElement(zeros(sim.nk+1), zeros(sim.nk+1)) for _ in 1:sim.nx+sim.numberOfGhostCells*2]
    temp_eulerianElements .= sim.eulerianElements
    for i in 1:sim.numberOfGhostCells
        sim.eulerianElements[i] = temp_eulerianElements[sim.nx + i]
        sim.eulerianElements[sim.nx + sim.numberOfGhostCells + i] = temp_eulerianElements[sim.numberOfGhostCells + i]
    end
end


# -----------------------------
# Update upstream vertex positions and compute corresponding modal coefficients.
# -----------------------------
function get_upstream_tn!(sim::SimulationState, ax::Function)
    # Step 1: Update vertex_star positions using Runge–Kutta.
    for i in 1:length(sim.vertices)
        new_coor = runge_kutta(sim.vertices[i].coor, sim.dt, ax)
        sim.vertex_star[i].coor = new_coor
        sim.vertex_star[i].id = ceil(Int, (new_coor - sim.xleft) / sim.dx)
    end

    # Step 2: Assemble Gauss–Lobatto nodes for upstream elements.
    for i in 1:sim.nx
        up = sim.upstreamElements[i]
        ee = sim.eulerianElements[i]
        # Set first and last upstream Gauss–Lobatto nodes from vertex_star.
        up.xgl_star[1] = sim.vertex_star[i].coor
        up.xgl_star[end] = sim.vertex_star[i+1].coor
        if sim.nk > 1
            # Set Eulerian element xgl nodes: first and last from vertices.
            ee.xgl[1] = sim.vertices[i].coor
            ee.xgl[end] = sim.vertices[i+1].coor
            # For interior nodes, use a placeholder linear mapping.
            for ii in 2:sim.nk
                ee.xgl[ii] = ee.xgl[1] + (0.5 + sim.gaussianQuadraturePoints[ii]) * (ee.xgl[end] - ee.xgl[1])
                up.xgl_star[ii] = runge_kutta(ee.xgl[ii], sim.dt, ax)
            end
        end

        # Step 3: Set upstream element endpoints and compute segments.
        up.point_origin = sim.vertex_star[i]
        up.point_end = sim.vertex_star[i+1]
        search_segment!(up, sim.xGrid, sim.numberOfGhostCells)

         # Step 4: Compute the modal coefficients via integration.
        
        out_coeffs = get_integral_pk!(up, sim)
        # println(out_coeffs)
        ee.umodal .= out_coeffs
    end
end

# -----------------------------
# Update the solution by copying upstream modal coefficients to Eulerian elements.
# -----------------------------
function update_solution!(sim::SimulationState)
    for i in 1:sim.nx
        sim.upstreamElements[i].umodal .= sim.eulerianElements[i].umodal
    end
end

# -----------------------------
# Compute error metrics by comparing the Eulerian solution to the exact solution.
# -----------------------------
function order_degree(sim::SimulationState, exact_function::Function)
    error1 = 0.0
    error2 = 0.0
    error3 = 0.0
    exact_ave = 0.0
    for i = 1 + sim.numberOfGhostCells:sim.nx + sim.numberOfGhostCells
        for ig = 1:sim.N
            xrg = sim.x[i] + (gausslegendre(sim.N))[1][ig]/2 * sim.dx
            a_1 = exact_function(xrg, sim.time_final)
            b_1 = ortho_poly1d(sim.eulerianElements[i].umodal, xrg, sim.x[i], sim.dx, sim.nk)
            exe =  a_1- b_1

            error1 += abs(exe)*gausslegendre(sim.N)[2][ig]/2
            error2 += exe^2*gausslegendre(sim.N)[2][ig]/2
            error3 = max(error3, abs(exe))
        end

    end
    error1 *= sim.dx
    error2 = sqrt(error2 * sim.dx)
    return error1, error2, error3
end

# -----------------------------
# --- Helper Functions ---
# Dummy implementation of search_segment! that fills point_inter, segment, and nsub.
function search_segment!(up::UpstreamElement, xGrid::Vector{Float64}, numberOfGhostCells::Int64)
    # Determine number of subintervals based on the indices.
    mx = up.point_end.id - up.point_origin.id
    up.point_inter = Vector{Vertex}(undef, mx + 2)
    up.point_inter[1] = up.point_origin
    up.point_inter[end] = up.point_end
    up.nsub = mx + 1
    if mx != 0
        for kk in 1:mx
            inter = up.point_origin.id + kk
            up.point_inter[kk+1] = Vertex(xGrid[inter + numberOfGhostCells], inter)
            # println(inter)
        end
    end
    # Create segments from the point_inter array.
    up.segment = Vector{Segment}(undef, mx + 1)
    for kk in 1:(mx + 1)
        up.segment[kk] = Segment(up.point_inter[kk], up.point_inter[kk+1], up.point_inter[kk].id)
        # println(up.point_inter[kk].id)
    end
end

# implementation of get_integral_pk! that computes the modal coefficients.
function get_integral_pk!(up::UpstreamElement, sim::SimulationState)
    # nk: polynomial degree; hence, there are nk+1 modes.
    nk = sim.nk  
    # Prepare the output array for the modal coefficients.
    out_umodal = zeros(nk+1)

    gau_t = zeros(nk+1)
    gau = gausslegendre(nk + 1)

    xi = zeros(nk+1)
    
    # Loop over test function modes: nm = 1,2,...,nk+1.
    for nm in 1:(nk+1)
        # Let "aa" represent the coefficients of the interpolated test function.
        aa = zeros(nk+1)
        if nm == 1
            # For the first mode: test function = 1.
            aa[1] = 1.0
            # (The remaining entries remain 0.)
        else
            if nm == 2
                # Compute the midpoint of the upstream Gauss–Lobatto nodes.
                xc_star = (up.xgl_star[1] + up.xgl_star[end]) / 2
                for igl in 1:(nk+1)
                    # Map the Gauss–Lobatto node into a reference coordinate.
                    xi[igl] = (up.xgl_star[igl] - xc_star) / (up.xgl_star[end] - up.xgl_star[1])
                end
                # Build the Vandermonde (or interpolation) matrix on the reference cell.
                amatrix = get_matrix(xi, nk)  # must be defined externally
                # Perform LU factorization of amatrix.
                store_L, store_U = doolittle(amatrix, nk+1)  # must be defined externally
            end
            # Evaluate the Legendre (or orthogonal) polynomials at the Gauss–Lobatto nodes.
            bb = zeros(nk+1)
            for kk in 1:(nk+1)
                bb[kk] = fle(nm-1, gausslobatto(nk+1)[1][kk]/2)  # fle must be defined; sim.gau_lob is assumed available.
            end
            # Solve the linear system using LU factors.
            aa = solve17(store_L, store_U, bb, nk+1)  # must be defined externally.
        end
        
        # Now compute the integral that yields the modal coefficient for this test function.
        sum_val = 0.0
        # Loop over each subinterval in the upstream element.
        for kk in 1:up.nsub
            # Get the index corresponding to the segment.            
            idx = up.segment[kk].id
            # Compute the quadrature nodes gau_t for this segment.
            for ig in 1:(nk+1)
                # For this segment, we assume an affine mapping:
                # x̄ = (porigin + pend)/2 + (pend - porigin)*node,
                porigin = up.segment[kk].porigin.coor
                pend    = up.segment[kk].pend.coor
                gau_t[ig] = (pend + porigin)/2 + (gau[1][ig]/2)*(pend - porigin)
            end
            
            # Determine the width of the entire upstream cell.
            dx_t = up.point_end.coor - up.point_origin.coor
            # Compute the midpoint of the upstream cell.
            xc_star = (up.point_origin.coor + up.point_end.coor) / 2
            st = 0.0
            # For each quadrature node in the test-function integration:
            for ig in 1:(nk+1)
                # Here we compute the product of two basis functions:
                # one coming from the reconstruction of the Eulerian solution in cell idx,
                # and the other from the test function (whose expansion is stored in aa).
                a_1 = ortho_poly1d(sim.eulerianElements[idx + sim.numberOfGhostCells].umodal, gau_t[ig], sim.x[idx + sim.numberOfGhostCells], sim.dx, nk)
                st += a_1 *
                      ortho_poly1d(aa, gau_t[ig], xc_star, dx_t, nk) *
                      (gau[2][ig]/2)
            end
            # Multiply by the segment length (scaled by the global dx).
            sum_val += st * (up.segment[kk].pend.coor - up.segment[kk].porigin.coor) / sim.dx

        end
        
        # Multiply by the appropriate correction coefficient (ai).
        out_umodal[nm] = sum_val * sim.coefficientsForCorrection[nm]
    end

    return out_umodal
end

function ortho_poly1d(a, x, xc, dx, k)
    if (k == 0)
        return a[1]
    elseif (k == 1)
        return a[1] + a[2] * (x - xc) / dx
    elseif (k == 2)
        return a[1] + a[2] * (x - xc) / dx + a[3] * (3*(x - xc)^2 - dx^2) / (2*dx^2)
    elseif (k == 3)
        return a[1] + a[2] * (x - xc) / dx + a[3] * (3*(x - xc)^2 - dx^2) / (2*dx^2) +
               a[4] * (5*(x - xc)^3 - 3*dx^2*(x - xc)) / (2*dx^3)
    end
    
end

function get_matrix(xi, nk)
    amatrix = zeros(nk+1, nk+1)
    if (nk == 1)
        amatrix[1,1] = 1.0
        amatrix[1,2] = xi[1]
        amatrix[2,1] = 1.0
        amatrix[2,2] = xi[2]
    elseif (nk == 2)
        amatrix[1,1] = 1.0
        amatrix[1,2] = xi[1]
        amatrix[1,3] = xi[1]^2 - 1/12
        amatrix[2,1] = 1.0
        amatrix[2,2] = xi[2]
        amatrix[2,3] = xi[2]^2 - 1/12
        amatrix[3,1] = 1.0
        amatrix[3,2] = xi[3]
        amatrix[3,3] = xi[3]^2 - 1/12 
    end
    return amatrix
end

function doolittle(A,N)
    L = zeros(N,N)
    U = zeros(N,N)
    for i in 1:N
        L[i,i] = 1.0
    end
    for i in 1:N
        for j in i:N
            sum = 0.0
            for k in 1:i-1
                sum += L[i,k]*U[k,j]
            end
            U[i,j] = A[i,j] - sum
        end
        for j in i+1:N
            sum = 0.0
            for k in 1:i-1
                sum += L[j,k]*U[k,i]
            end
            L[j,i] = (A[j,i] - sum) / U[i,i]
        end
    end
    return L, U
    
end

function solve17(L, U, b, N)
    y = zeros(N)
    x = zeros(N)
    # Solve Ly = b
    y[1] = b[1]/L[1,1]
    for i = 2:N
        y[i] = b[i]
        for j = 1:i-1
            y[i] -= L[i,j]*y[j]
        end
        y[i] /= L[i,i]
    end
    # Solve Ux = y
    x[N] = y[N]/U[N,N]
    for i = N-1:-1:1
        x[i] = y[i]
        for j = i+1:N
            x[i] = x[i] - U[i,j]*x[j]
        end
        x[i] /= U[i,i]
    end
    return x
end

function fle(k,x)
    if(k == 0)
        return 1.0
    elseif(k == 1)
        return x
    elseif(k == 2)
        return x^2 - 1/12
    elseif(k == 3)
        return x^3 - 3/20*x
    elseif(k == 4)
        return x^4 - 6/35*x^2 + 3/140
    end
end


# -----------------------------
# Master function: sldg1d
# -----------------------------
function sldg1d(params::SimulationParameters, ax::Function, fun_init::Function, exact_function::Function)
    # Create a SimulationState.
    sim = SLDG1d.SimulationState(
        params.nx,                        # nx
        params.nk,                        # nk
        params.numberOfGhostCells,        # numberOfGhostCells
        params.N,                         # N (quadrature points)
        0.0,                              # dx (to be set in setup!)
        0.0,                              # pi
        0.0,                              # eps
        zeros(params.nk+1),               # coefficientsForCorrection
        params.xleft,                     # xleft
        params.xright,                    # xright
        0.0,                              # dt
        0.0,                              # time
        params.time_final,                # time_final
        params.cfl,                       # cfl
        Float64[],                        # x (grid cell centers, allocated in setup!)
        Float64[],                        # xGrid (grid including ghost cells)
        Vertex[],                         # vertices (cell boundaries)
        Vertex[],                         # vertex_star (upstream cell boundaries)
        EulerianElement[],                # eulerianElements (one per cell)
        UpstreamElement[]                 # upstreamElements (one per cell)
    )
    
    # Initialize simulation constants and quadrature arrays.
    parameters!(sim)
    
    # Set up grid, vertices, and element arrays.
    setup!(sim, params)
    
    # Initialize the Eulerian elements with the initial condition.
    initialize!(sim, params, fun_init)
    
    # Main time-stepping loop.
    while sim.time < sim.time_final
        update_dt(sim)
        boundary!(sim)
        get_upstream_tn!(sim, ax)
        update_solution!(sim)
    end
    
    
    # Compute error metrics.
    er1, er2, er3 = order_degree(sim, exact_function)
    
    # Package outputs.
    # Create a matrix (nx x (nk+1)) to hold the modal coefficients from Eulerian elements.
    solution = zeros(sim.nx, sim.nk+1)
    for i in 1:sim.nx
        solution[i, :] = sim.eulerianElements[i].umodal
    end
    
    # Compute cell centers from the vertices.
    grid = [0.5 * (sim.vertices[i].coor + sim.vertices[i+1].coor) for i in 1:sim.nx]
    
    return solution, grid, er1, er2, er3
end


end

