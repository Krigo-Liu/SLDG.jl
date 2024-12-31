# Set the Golable variable

## Set the whole domain
xLeft = 0.0
xRight = 2 .* pi

## Set the number of ghost elements
numberGhostElements = 11

## Set the time relative parameters
CFL_number = 0.1
timeFinal = 20.0

## Set the p\^th DG solution space
p = 1

# Set the structure of the problem

# Set the initial condition of PDEs
function sinInitialCondition(real :: x)
    return sin(x)
end

# Set the flux of PDEs???
function flux_X(x::Float64)
    return 1.0
    # return sin(x)
end

# Set the exact solution of PDEs
function exact_Sol(x, t)
    return sin(x - t)
    # return sin( 2.0 * atan(exp(-t) * tan(x / 2.0))) / sin(x)
end

for kkkk in 1:7

    ## The number of elements in the mesh
    nx = 10*2^kkkk
    dx = (xLeft - xRight) / nx

    for i = 1 - numberGhostElements:nx + 1 + numberGhostElements in 
        x[i] = xLeft + (i - 0.5) * dx
    end

    for i = 1 - numberGhostElements:nx + 1 + numberGhostElements in 
        xGrid[i] = xLeft + (i - 1.0) * dx
    end

    for i = 1:nx+1 in 
        
    end

    
end








