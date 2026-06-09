using QuadGK

n = 11 
x,w,wg = kronrod(n,-1,1)

function trace_finite(a,b)
    # parametrisation of finite straight line from a to b
    u -> 0.5*((b+a) + (b-a)*u) # :: Function
end

# Default parameters for Gauss-Kronrod 
# Nodes and weights for G7/K15 precomputed for (-1,1) using QuadGK.jl
const x_gk  = [-0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993945, -0.5860872354676911, -0.4058451513773972, -0.2077849550078985, 0.0, 0.2077849550078985, 0.4058451513773971, 0.5860872354676911, 0.7415311855993945, 0.8648644233597691, 0.9491079123427584, 0.9914553711208125]
const w_gk  = [0.022935322010529256, 0.06309209262997842, 0.10479001032225017, 0.14065325971552592, 0.16900472663926788, 0.19035057806478559, 0.20443294007529877, 0.20948214108472793, 0.20443294007529877, 0.19035057806478559, 0.16900472663926788, 0.14065325971552592, 0.10479001032225017, 0.06309209262997842, 0.022935322010529256]
const wg_gk = [0.12948496616886981, 0.2797053914892767, 0.38183005050511887, 0.41795918367346907, 0.38183005050511887, 0.2797053914892767, 0.12948496616886981]

function eval_gk(f, h :: Function, dh)
    fx = f.(h.(x_gk))
    I  = dh * sum(w_gk .* fx) 
    Ig = dh * sum(wg_gk .* fx[2:2:end])
    aerr = abs(I-Ig) # absolute error
    rerr = aerr / abs(I) # relative error
    return I, aerr, rerr
end

function eval_gk(f,a,b)
    h  = trace_finite(a,b)
    dh = 0.5 * (b-a)
    eval_gk(f,h,dh)
end

function do_quadgk(f,a,b; rtol = 1e-8, atol = 1e-8)
    # @show (a,b)
    val, abserror, relerror = eval_gk(f,a,b)
    if relerror > rtol && abserror > atol
        mid = 0.5*(a+b)
        # should return the sum of the two intervals
        val1, e1 = do_quadgk(f,a,mid) 
        val2, e2 = do_quadgk(f,mid,b)
        return val1+val2, e1+e2
    end
    return val, relerror
end

f(x) = x^2
f(x) = cos(100*x)
f(x) = 1/sqrt(x)

a = 0
b = 3
rtol = 1e-11
@time exact = quadgk(f,a,b, rtol)
@time mygk = do_quadgk(f,a,b)

mygk[1]-exact[1]

@time do_quadgk(f,a,b)
@profview do_quadgk(f,a,b)