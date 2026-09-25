module F03VectorCalculus

export scalar_field, vector_field, gradient_scalar, curl_vector, laplacian_scalar,
    divergence_vector, gradient_divergence_vector, laplacian_vector

scalar_field(point) = sin(point[1]) * cos(point[2]) * exp(point[3])

vector_field(point) = (
    sin(point[2]) * exp(2point[3]) + exp(point[1]),
    sin(point[3]) * exp(2point[1]) + exp(point[2]),
    sin(point[1]) * exp(2point[2]) + exp(point[3]),
)

function validate_point(point)
    point isa Tuple && length(point) == 3 ||
        throw(ArgumentError("点は3要素のタプルで指定してください"))
    all(value -> value isa Real && isfinite(value), point) ||
        throw(ArgumentError("点の座標は有限な実数にしてください"))
    nothing
end

function gradient_scalar(point)
    validate_point(point)
    (
        cos(point[1]) * cos(point[2]) * exp(point[3]),
        -sin(point[1]) * sin(point[2]) * exp(point[3]),
        sin(point[1]) * cos(point[2]) * exp(point[3]),
    )
end

function curl_vector(point)
    validate_point(point)
    (
        2sin(point[1]) * exp(2point[2]) - cos(point[3]) * exp(2point[1]),
        2sin(point[2]) * exp(2point[3]) - cos(point[1]) * exp(2point[2]),
        2sin(point[3]) * exp(2point[1]) - cos(point[2]) * exp(2point[3]),
    )
end

function laplacian_scalar(point)
    validate_point(point)
    -scalar_field(point)
end

function divergence_vector(point)
    validate_point(point)
    sum(exp, point)
end

function gradient_divergence_vector(point)
    validate_point(point)
    exp.(point)
end

function laplacian_vector(point)
    validate_point(point)
    x, y, z = point
    (3sin(y) * exp(2z) + exp(x),
     3sin(z) * exp(2x) + exp(y),
     3sin(x) * exp(2y) + exp(z))
end

end
