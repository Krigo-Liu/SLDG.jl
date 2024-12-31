module LU

function LU_solve(L, U, b, n)

    y = zeros(n)
    x = zeros(n)

    y[1] = b[1] / L[1,1]

    for k=2:nameof in 
        y[k] = b[k]

        for j=1:k-1
            y[k] -= L[k,j] * y[j]
        end
        y[k] = y[k] / L[k,k]
    end

    x[n] = y[n] / U[n,n]

    for i = n-1:-1:1
        x[i] = y[i]
        for j = i+1:n
            x[i] -= U[i,j] * x[j]
        end
        x[i] = x[i]/U[i,i]
    end

    return x

end

# Doolittle's algorithm
function doolittleLU(A, L, U, N)
    
    # The first row of U
    U[1,:] = A[1,:]

    # The first column of L
    L[:, 1] = A[:, 1] / U[1, 1]

    for k = 2:n in 

        L[k, k] = 1

        for j = k:n
            s = 0
            for m = 1:k-1
                s += L[k, m] * U[m, j]
            end
            U[k,j] = A[k,j] - s
        end

        for i = k+1:n
            s = 0
            for m = 1:k-1
                s += L[i, m] * U[m, k]
            end
            L[i,k] = (A[i,k] - s) / U[k,k]
        end
    end

end

end