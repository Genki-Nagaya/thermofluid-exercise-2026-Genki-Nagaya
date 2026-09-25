if !isdefined(Main, :F03VectorCalculus)
    include(joinpath(@__DIR__, "F03.jl"))
end

module F04NumericalDifferentiation

using Main.F03VectorCalculus: scalar_field, vector_field, gradient_scalar, curl_vector,
    divergence_vector, gradient_divergence_vector, laplacian_vector

export forward_difference,
    backward_difference,
    centered_difference,
    convergence_study,
    centered_partial,
    curl_gradient_residual,
    divergence_curl_residual,
    product_divergence_residual,
    curl_curl_residual,
    verify_vector_identities

function validate_scalar_input(x, h)
    x isa Real && !(x isa Bool) && isfinite(x) ||
        throw(ArgumentError("xは有限な実数にしてください"))
    h isa Real && !(h isa Bool) && isfinite(h) && h > 0 ||
        throw(ArgumentError("hは有限な正の実数にしてください"))
    nothing
end

function validate_point(point)
    point isa Tuple && length(point) == 3 ||
        throw(ArgumentError("点は3要素のタプルで指定してください"))
    all(value -> value isa Real && !(value isa Bool) && isfinite(value), point) ||
        throw(ArgumentError("点の座標は有限な実数にしてください"))
    nothing
end

function forward_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 前進差分商を実装する。
    zero(float(x + h))
end

function backward_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 後退差分商を実装する。
    zero(float(x + h))
end

function centered_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 中心差分商を実装する。
    zero(float(x + h))
end

function convergence_study(f, derivative, x, spacings)
    x isa Real && !(x isa Bool) && isfinite(x) ||
        throw(ArgumentError("xは有限な実数にしてください"))
    spacings isa AbstractVector && length(spacings) >= 2 ||
        throw(ArgumentError("spacingsには少なくとも2個の値が必要です"))
    all(h -> h isa Real && !(h isa Bool) && isfinite(h) && h > 0, spacings) ||
        throw(ArgumentError("spacingsは有限な正の実数にしてください"))
    all(index -> spacings[index] > spacings[index + 1], 1:(length(spacings) - 1)) ||
        throw(ArgumentError("spacingsは狭義単調減少にしてください"))

    exact = derivative(x)
    forward_errors = [abs(forward_difference(f, x, h) - exact) for h in spacings]
    backward_errors = [abs(backward_difference(f, x, h) - exact) for h in spacings]
    centered_errors = [abs(centered_difference(f, x, h) - exact) for h in spacings]
    ratios(errors) = errors[1:(end - 1)] ./ errors[2:end]
    forward_ratios = ratios(forward_errors)
    backward_ratios = ratios(backward_errors)
    centered_ratios = ratios(centered_errors)
    (;
        spacings,
        forward_errors,
        backward_errors,
        centered_errors,
        forward_ratios,
        backward_ratios,
        centered_ratios,
        forward_orders = log2.(forward_ratios),
        backward_orders = log2.(backward_ratios),
        centered_orders = log2.(centered_ratios),
    )
end

function centered_partial(f, point, axis, h; differentiator = centered_difference)
    validate_point(point)
    h isa Real && !(h isa Bool) && isfinite(h) && h > 0 ||
        throw(ArgumentError("hは有限な正の実数にしてください"))
    axis isa Integer && !(axis isa Bool) && axis in 1:3 ||
        throw(ArgumentError("axisは1、2、3のいずれかにしてください"))
    slice(value) = f(ntuple(index -> index == axis ? value : point[index], 3))
    differentiator(slice, point[axis], h)
end

# The supplied analytic gradient is differentiated externally with the student's centered difference.
# The inner gradient is not student work in this residual calculation.
function curl_gradient_residual(point, h)
    validate_point(point)
    (
        centered_partial(p -> gradient_scalar(p)[3], point, 2, h) -
        centered_partial(p -> gradient_scalar(p)[2], point, 3, h),
        centered_partial(p -> gradient_scalar(p)[1], point, 3, h) -
        centered_partial(p -> gradient_scalar(p)[3], point, 1, h),
        centered_partial(p -> gradient_scalar(p)[2], point, 1, h) -
        centered_partial(p -> gradient_scalar(p)[1], point, 2, h),
    )
end

# The supplied analytic curl is differentiated externally with the student's centered difference.
# The inner curl is not student work in this residual calculation.
function divergence_curl_residual(point, h)
    validate_point(point)
    sum(centered_partial(p -> curl_vector(p)[axis], point, axis, h) for axis in 1:3)
end

# Differentiate the product field numerically and evaluate the right-hand side analytically.
function product_divergence_residual(point, h)
    validate_point(point)
    lhs = sum(centered_partial(p -> scalar_field(p) * vector_field(p)[axis],
        point, axis, h) for axis in 1:3)
    rhs = sum(vector_field(point) .* gradient_scalar(point)) +
        scalar_field(point) * divergence_vector(point)
    lhs - rhs
end

function curl_curl_residual(point, h)
    validate_point(point)
    lhs = (
        centered_partial(p -> curl_vector(p)[3], point, 2, h) -
        centered_partial(p -> curl_vector(p)[2], point, 3, h),
        centered_partial(p -> curl_vector(p)[1], point, 3, h) -
        centered_partial(p -> curl_vector(p)[3], point, 1, h),
        centered_partial(p -> curl_vector(p)[2], point, 1, h) -
        centered_partial(p -> curl_vector(p)[1], point, 2, h),
    )
    lhs .- (gradient_divergence_vector(point) .- laplacian_vector(point))
end

function verify_vector_identities(n)
    n isa Integer && !(n isa Bool) && n >= 5 && isodd(n) ||
        throw(ArgumentError("nは5以上の奇数にしてください"))
    coordinates = range(-1.0, 1.0; length = Int(n))
    h = step(coordinates)
    curl_gradient = 0.0
    divergence_curl = 0.0
    product_divergence = 0.0
    curl_curl = 0.0
    # Compare maxima on a common physical region as the grid is refined.
    interior = filter(x -> abs(x) <= 0.75, coordinates[2:(end - 1)])
    for x in interior, y in interior, z in interior

        point = (x, y, z)
        curl_gradient = max(curl_gradient, maximum(abs, curl_gradient_residual(point, h)))
        divergence_curl = max(divergence_curl, abs(divergence_curl_residual(point, h)))
        product_divergence = max(product_divergence, abs(product_divergence_residual(point, h)))
        curl_curl = max(curl_curl, maximum(abs, curl_curl_residual(point, h)))
    end
    (; curl_gradient, divergence_curl, product_divergence, curl_curl)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    reference_function(x) = sin(x) * exp(x)
    reference_derivative(x) = exp(x) * (sin(x) + cos(x))
    spacings = [0.2, 0.1, 0.05, 0.025]
    println(
        F04NumericalDifferentiation.convergence_study(
            reference_function,
            reference_derivative,
            0.4,
            spacings,
        ),
    )
    println("n=9: ", F04NumericalDifferentiation.verify_vector_identities(9))
    println("n=17: ", F04NumericalDifferentiation.verify_vector_identities(17))
end
