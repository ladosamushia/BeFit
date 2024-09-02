using DelimitedFiles
using Interpolations

function PL(k, PLint)
    return PLint(k)
end

function PLint(k, Pk)
    return linear_interpolation(k, Pk)
end

function PLint_from_file(filename)
    Pkk = readdlm(filename, comments=true)
    k = append!([0.0], Pkk[:,1])
    Pk = append!([0.0], Pkk[:,2])
    return PLint(k, Pk)
end