# Function for Legendre polynomials

function legendrePoly(n, x)
    if n == 0
        return 1.0
    
    elseif n == 1
        return x
    
    elseif n == 2
        return x*x - 1.0/12.0
    
    elseif n == 3
        return x*x*x - 0.15 * x
    
    elseif n == 4
        return (x^2 - 3.0 / 14.0) *x*x + 3.0/560.0
    end
end