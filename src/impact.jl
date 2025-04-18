"""
"""
function augmented_impact_map(
    manip::ManipulatorEquation,
    e::Float64,
    q::Vector{Float64}
)::Matrix{Float64}
    Minv, Jinv, Λinv = manipulator_inverses(manip, x)
    M = manip.M(q)
    J = manip.J(q)
    return J' * (e*Λinv*J - Jinv*M)
end

"""
"""
function impact_law(
    M::Matrix{Float64},
    J::Matrix{Float64},
    v::Vector{Float64}
)::Vector{Float64}
    # Assume that all constraint terms are associated with impact
    Minv = M \ I
    return (I + Minv*J')*v
end
