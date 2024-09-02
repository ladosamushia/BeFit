include("RedshiftSpaceKernels.jl")

using HCubature

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

function B222(k1, k2, f, bias, PL, kmax)
    fint = x -> B222int(k1, k2, x, f, bias, PL)
    result = hcubature(fint, -kmax, kmax, rtol=1e-12, maxevals=1000000)
    return result
end

function B321I()
end

function B321II()

end

function B411()

end
