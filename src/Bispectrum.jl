include("RedshiftSpaceKernels.jl")

using HCubature
using Combinatorics
using Cuba

function B211(k1, k2, f, bias, PL)
    result = 2*Z2(k1, k2, f, bias)*Z1(k1, f, bias)*Z1(k2, f, bias)*PL(qlen(k1))*PL(qlen(k2))
    result += 2*Z2(k2, k3, f, bias)*Z1(k2, f, bias)*Z1(k3, f, bias)*PL(qlen(k2))*PL(qlen(k3))
    result += 2*Z2(k3, k1, f, bias)*Z1(k3, f, bias)*Z1(k1, f, bias)*PL(qlen(k3))*PL(qlen(k1))
    return result
end

function B222int(k1, k2, q, f, bias, PL)
    result = 8*Z2(k1 + q, -q, f, bias)*Z2(k1 + q, k2 - q, f, bias)*Z2(k2 - q, q, f, bias)*PL(qlen(q))*PL(qlen(k1 + q))*PL(qlen(k2 - q))
    return result
end

function B222intSph(k1, k2, qsp, f, bias, PL)
    q = [qsp[1]*sin(qsp[2])*cos(qsp[3]), qsp[1]*sin(qsp[2])*sin(qsp[3]), qsp[1]*cos(qsp[2])]
    result = 8*Z2(k1 + q, -q, f, bias)*Z2(k1 + q, k2 - q, f, bias)*Z2(k2 - q, q, f, bias)*PL(qlen(q))*PL(qlen(k1 + q))*PL(qlen(k2 - q))
    return result
end

function B222intAng(k1, k2, q, omega, f, bias, PL)
    
    q = [q*sin(omega[1])*cos(omega[2]), q*sin(omega[1])*sin(omega[2]), q*cos(omega[1])]
    result = 8*Z2(k1 + q, -q, f, bias)*Z2(k1 + q, k2 - q, f, bias)*Z2(k2 - q, q, f, bias)*PL(qlen(q))*PL(qlen(k1 + q))*PL(qlen(k2 - q))
    return result
end

function B222(k1, k2, f, bias, PL, kmax)
    result = (0, 0)
    for pair in combinations((k1,k2,-k1-k2), 2)
        fint = x -> B222int(pair[1], pair[2], x, f, bias, PL)
        result = result .+ hcubature(fint, -kmax, kmax)
    end
    return result
end

function B222sp(k1, k2, f, bias, PL, kmax)
    result = (0, 0)
    for pair in combinations((k1,k2,-k1-k2), 2)
        fint = (x, func) -> func[1] = B222intSph(pair[1], pair[2], x./[kmax, pi, 2*pi], f, bias, PL)*(x[1]/kmax)^2*sin(x[2]/pi)
        result = result .+ suave(fint, 3, rtol=1e-3).integral[1]
    end
    return result
end

function B222ang(k1, k2, q, f, bias, PL)
    result = 0
    for pair in combinations((k1,k2,-k1-k2), 2)
        fint = (x, func) -> func[1] = B222intAng(pair[1], pair[2], q, x./[pi, 2*pi], f, bias, PL)*sin(x[1]/pi)
        result = result .+ suave(fint, 2, rtol=1e-3).integral[1]
    end
    return result
end
function B321I()
end

function B321II()

end

function B411()

end
