# k1 and k2 need to be 3-element Vectors
# k1 = [1, 2, 3]

function alpha(k1, k2)
    return (k1 + k2)'*k1/(k1'*k1)
end

function beta(k1, k2)
    return (k1 + k2)'*(k1 + k2)*k1'*k2/2/(k1'*k1)/(k2'*k2)
end
# n must be an integer. k must be a matrix of size 3xn
# k = [k1 k2 k3 ...]
function Fkernel(n, k)
    F = 0
    if n == 1
        F = 1
    else
        for m in 1:n-1
            k1 = vec(sum(k[:,1:m], dims=2))
            k2 = vec(sum(k[:,m+1:end], dims=2))
            F += Gkernel(m, k[:,1:m])/(2*n + 3)/(n - 1)*((2*n + 1)*alpha(k1, k2)*Fkernel(n - m, k[:,m+1:end]) + 2*beta(k1, k2)*Gkernel(n - m, k[:,m+1:end]))
        end
    end
    return F
end

function Gkernel(n, k)
    G = 0
    if n == 1
        G = 1
    else
        for m in 1:n-1
            k1 = vec(sum(k[:,1:m], dims=2))
            k2 = vec(sum(k[:,m+1:end], dims=2))
            G += Gkernel(m, k[:,1:m])/(2*n + 3)/(n - 1)*(3*alpha(k1, k2)*Fkernel(n - m, k[:,m+1:end]) + 2*n*beta(k1, k2)*Gkernel(n - m, k[:,m+1:end]))
        end
    end
    return G
end
