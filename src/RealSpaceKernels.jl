include("EFTkernels.jl")

function mu(q1, q2)
    return q1'*q2/sqrt(q1'*q1)/sqrt(q2'*q2)
end

function kappa(q1, q2)
    return mu(q1, q2)^2 - 1
end

function L(q1, q2, q3)
    return 2*mu(q1, q2)*mu(q2, q3)*mu(q3, q1) - mu(q1, q2)^2 - mu(q2, q3)^2 - mu(q3, q1)^2 + 1
end

function M(q1, q2, q3, q4)
    return mu(q1, q2)*mu(q2, q3)*mu(q3, q4)*mu(q4, q1)
end

function K1(q1, bias)
    result = bias["b1"]
    return result
end

function K2(q1, q2, bias)
    result = bias["b1"]*Fkernel(2, [q1 q2]) 
    result += bias["b2"]/2 + bias["gamma2"]*kappa(q1, q1 + q2)
    return result
end

function K3(q1, q2, q3, bias)
    result = bias["b1"]*Fkernel(3, [q1 q2 q3]) 
    result += bias["b2"]*Fkernel(2, [q1 q2]) + 2*bias["gamma2"]*kappa(q1, q2)*Gkernel(2, [q2 q3])
    result += bias["b3"]/6 + bias["gamma2x"]*kappa(q1, q2) + bias["gamma3"]*L(q1, q2, q3) + bias["gamma21"]*kappa(q1, q2)*Gkernel(3, [q1 q2 q3])
    return result
end

function K4(q1, q2, q3, q4, bias)
    result = bias["b1"]*Fkernel(4, [q1 q2 q3 q4])
    result += bias["b2"]/2*(Fkernel(2, [q1 q2])*Fkernel(2, [q3 q4]) + 2*Fkernel(3, [q1 q2 q3])) + bias["gamma2"]*(kappa(q1 + q2, q3 + q4)*Gkernel(2, [q1 q2])*Gkernel(2, [q3 q4]) + 2*kappa(q1 + q2 + q3, q4)*Gkernel(3, [q1 q2 q3]))
    result += bias["b3"]/2*Fkernel(2, [q1 q2]) + bias["gamma2x"]*(2*kappa(q1 + q2, q3)*Gkernel(2, [q1 q2]) + kappa(q3, q4)*Fkernel(2, [q1 q2])) + 3*bias["gamma3"]*L(q1, q2, q3 + q4)*Gkernel(2, [q3 q4]) + bias["gamma21"]*(kappa(q1 + q2, q3 + q4)*kappa(q1, q2)*Fkernel(2, [q3 q4]) + 2*kappa(q1 + q2 + q3, q4)*kappa(q1 + q2, q3)*Fkernel(2, [q1 q2]))
    result += bias["gamma21x"]*kappa(q1, q2 + q3)*kappa(q2, q3) + bias["gamma211"]*L(q1, q2, q3 + q4)*kappa(q3, q4) + bias["gamma22"]*kappa(q1 + q2, q3 + q4)*kappa(q1, q2)*kappa(q3, q4) + bias["gamma31"]*(1/18*kappa(q1, q2 + q3 + q4)*(15/7*kappa(q2 + q3, q4)*kappa(q2, q3) - L(q2, q3, q4)) + 1/14*(M(q1, q2 + q3, q4, q2 + q3 + q4) - M(q1, q2 + q3 + q4, q2 + q3, q4))*kappa(q2, q3))
    return result 
end
