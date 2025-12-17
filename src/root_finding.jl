
function roots_trig_polynomial(a,b)

    """ 
    We are following "Computing the zeros, maxima and inflection points
        of Chebyshev, Legendre and Fourier series: solving
        transcendental equations by spectral interpolation
        and polynomial rootfinding" (Boyd, 2007)
    """

    # assume a vector of coefficients of cos(nθ) and b vector of coefficients of sin(nθ)
    # length(a) = 2N + 1
    @assert length(a) == length(b) "Coefficient vectors must be of the same length"
    # @assert abs(b[1]) ≈ 0.0 "The first coefficient of b_k (sin(kθ)) must be zero"
    J = length(a) - 1

    # compute Foorier companion matrix  
    h = zeros(ComplexF64, 2J+1)
    h[1:J]      = [a[J-j+1] + im*b[J-j+1] for j = 0:J-1]
    h[J+1]      = 2*a[1]
    h[J+2:2J+1] = [a[j-J+1] - im*b[j-J+1] for j = J+1:2J]

    B = zeros(ComplexF64, 2J, 2J)
    B[1:end-1, 2:end] = I(2J-1)
    B[end, :] = - h[1:end-1] / (a[end] - im*b[end]) 

    # convert eigenvalues to roots
    z = eigvals(B)
    θ = angle.(z) - im*log.(abs.(z)) 
    return θ
end