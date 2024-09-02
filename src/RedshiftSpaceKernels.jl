include("RealSpaceKernels.jl")

function mu(q)
    return q[3]/sqrt(q'*q)
end

function qlen(q)
    return sqrt(q'*q)
end

function Z1(q1, f, bias)
    result = K1(q1, bias) + f*mu(q1)^2
    return result 
end

function Z2(q1, q2, f, bias)
    result = K2(q1, q2, bias) + f*mu(q1 + q2)^2*Gkernel(2, [q1 q2]) 
    result += f*mu(q1 + q2)*qlen(q1 + q2)/2*K1(q1, bias)*(mu(q1)/qlen(q1) + mu(q2)/qlen(q2)) 
    result += (f*mu(q1 + q2)*qlen(q1 + q2))^2/2*mu(q1)/qlen(q1)*mu(q2)/qlen(q2)
    return result 
end

function Z3(q1, q2, q3, f, bias)
    result = K3(q1, q2, q3, bias) + f*mu(q1 + q2 + q3)^2*Gkernel(3, [q1 q2 q3])
    result += f*mu(q1 + q2 + q3)*qlen(q1 + q2 + q3)*(mu(q1 + q2)/(q1 + q2)*K1(q1, bias)*G2(q1, q2) + mu(q3)/qlen(q3)*K2(q1, q2, bias))
    result += (f*mu(q1 + q2 + q3)*qlen(q1 + q2 + q3))^2/2*(2*mu(q1 + q2)/qlen(q1 + q2)*mu(q3)/qlen(q3)*G2(q1, q2, bias) + mu(q1)/qlen(q1)*mu(q2)/qlen(q2)*K1(q1, bias))
    result += (f*mu(q1 + q2 + q3)*qlen(q1 + q2 + q3))^3/6*mu(q1)/qlen(q1)*mu(q2)/qlen(q2)*mu(q3)/qlen(q3)
    return result 
end

function Z4(q1, q2, q3, q4, f, bias)
    result = K4(q1, q2, q3, q4, bias) + f*mu(q1 + q2 + q3 + q4)^2*Gkernel(4, [q1 q2 q3 q4])
    result += f*mu(q1 + q2 + q3 + q4)*qlen(q1 + q2 + q3 + q4)*(mu(q1 + q2 + q3)/qlen(q1 + q2 + q3)*Gkernel(3,[q1 q2 q3]) + mu(q4)/qlen(q4)*K3(q1, q2, q3, bias) + mu(q1 + q2)/qlen(q1 + q2)*Gkernel(2, [q1 q2])*K2(q1, q2, bias))
    result += (f*mu(q1 + q2 + q3 + q4)*qlen(q1 + q2 + q3 + q4))^2/2*(2*mu(q1 + q2 + q3)/qlen(q1 + q2 + q3)*mu(q4)/qlen(q4)*Gkernel(3, [q1 q2 q3]) + mu(q1 + q2)/qlen(q1 + q2)*mu(q3 + q4)/qlen(q3 + q4)*Gkernel(2, [q1 q2])*Gkernel(2, [q3 q4]) + 2*mu(q1 + q2)/qlen(q1 + q2)*mu(q3)/qlen(q3)*Gkernel(2, [q1 q2]) + mu(q1)/qlen(q1)*mu(q2)/qlen(q2)*Kkernel(2, [q3 q4]))
    result += (f*mu(q1 + q2 + q3 + q4)*qlen(q1 + q2 + q3 + q4))^3/6*(3*mu(q1 + q2)/qlen(q1 + q2)*mu(q3)/qlen(q3)*mu(q4)/qlen(q4)*Gkernel(2, [q1 q2]) + mu(q1)/qlen(q1)*mu(q2)/qlen(q2)*mu(q3)/qlen(q3))
    result += (f*mu(q1 + q2 + q3 + q4)*qlen(q1 + q2 + q3 + q4))^4/24*mu(q1)/qlen(q1)*mu(q2)/qlen(q2)*mu(q3)/qlen(q3)*mu(q4)/qlen(q4)
    return result 
end