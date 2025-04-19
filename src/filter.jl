"""
"""
mutable struct SafetyFilter
    nu::Int
    ϵ::Real
    Φ::Function
    Φ̇ub::Function
    flow::ControlAffineFlow
    #model::OSQP.Model
    #settings::NamedTuple
end

#= function SafetyFilter(
    nu::Int,
    ϵ::Real,
    Φ::Function,
    Φ̇ub::Function,
    flow::ControlAffineFlow
)::SafetyFilter
    # Get safety index dimension
    return SafetyFilter(nu, ϵ, Φ, Φ̇ub, flow)

    #model = OSQP.Model()
    #P = Matrix{<:Real}(2*I(nu))
    #A = zeros(nΦ, nu)
    #q = zeros(nu)
    #l = zeros(nu)
    #OSQP.setup!(model; P=P, A=A, q=q, l=l, u=l, settings...)
    #return SafetyFilter(nΦ, Φ, Φ̇ub, flow, model, settings)
end =#


    ## Compute safe control set bounds
    #δΦ = ForwardDiff.jacobian(filter.Φ, x)
    #A = δΦ*filter.flow.actuated(x)
    #q = -2*u
    #lb = -Inf*ones(filter.nΦ)
    #ub = filter.Φ̇ub(x) - δΦ*filter.flow.unactuated(x)

    ## Solve QP
    #OSQP.update!(filter.model, Ax=A; q=q, l=lb, u=ub)
    #OSQP.update_settings!(filter.model; filter.settings)
    #results = OSQP.solve!(filter.model)
    #return results.x


"""
"""
function (filter::SafetyFilter)(
    Φargs::Tuple,
    Φ̇args::Tuple,
    x::Vector{<:Real},
    u::Vector{<:Real};
    verbose::Bool = false
)::Vector{<:Real}
    # Quadratic program
    uvar = Variable(filter.nu)
    prob = minimize(square(uvar - u))

    # Compute safe control set bounds
    Φ̇ub = filter.Φ̇ub(Φ̇args..., x)
    δΦ = ForwardDiff.jacobian(δx -> filter.Φ(Φargs..., δx), x)
    f = filter.flow.actuated(x)
    g = filter.flow.unactuated(x)
    prob.constraints = [δΦ*(f + g*uvar) <= Φ̇ub .+ filter.ϵ]

    # Solve QP/SOCP
    solve!(prob, Clarabel.Optimizer; silent = !verbose)
    result = evaluate(uvar)
    if typeof(result) <: Real
        return [result]
    else
        return result
    end
end
