"""
"""
mutable struct SafetyFilter
    nu::Int
    ϵ::Real
    Φ::Function
    Φ̇ub::Function
    flow::ControlAffineFlow
    uvar::Vector{VariableRef}
    model::Model
    #settings::NamedTuple
end

function SafetyFilter(
    nu::Int,
    ϵ::Real,
    Φ::Function,
    Φ̇ub::Function,
    flow::ControlAffineFlow
)::SafetyFilter
    model = Model(Clarabel.Optimizer)
    @variable(model, uvar[1:nu])
    return SafetyFilter(nu, ϵ, Φ, Φ̇ub, flow, uvar, model)
end

"""
"""
function (filt::SafetyFilter)(
    Φargs::Tuple,
    Φ̇args::Tuple,
    x::Vector{<:Real},
    uref::Vector{<:Real};
    verbose::Bool = false
)::Vector{<:Real}
    # Activate filter if safety constraint is violated
    #=use_filter = false
    Φ = filt.Φ(Φargs..., x)
    for val in Φ
        if val > 0
            use_filter = true
            break
        end
    end
    use_filter = true

    if !use_filter
        return uref
    else=#
    # Compute safe control set bounds
    Φ̇ub = filt.Φ̇ub(Φ̇args..., x)
    δΦ = ForwardDiff.jacobian(δx -> filt.Φ(Φargs..., δx), x)
    f = reshape(filt.flow.unactuated(x), length(x))
    g = reshape(filt.flow.actuated(x), length(x), filt.nu)
    cexpr = δΦ*(f + g*filt.uvar)

    # Set objective
    obj = filt.uvar'*filt.uvar - 2*uref'*filt.uvar + uref'*uref
    @objective(filt.model, Min, obj)

    # Set constraints
    for (i, ub) in enumerate(Φ̇ub)
        if !isinf(ub)
            @constraint(filt.model, cexpr[i] <= ub + filt.ϵ)
        end
    end

    # Set solver options
    if verbose
        set_optimizer_attribute(filt.model, "verbose", true)
    else
        set_optimizer_attribute(filt.model, "verbose", false)
    end

    # Solve
    optimize!(filt.model)
    u = value.(filt.uvar)

    # Clear constraints
    for c in all_constraints(filt.model; include_variable_in_set_constraints = false)
        delete(filt.model, c)
    end

    # Return solution as a vector
    if typeof(u) <: Real
        return [u]
    else
        return u
    end
end


"""
"""
mutable struct ImpactAwareSafetyFilter
    q::Real
    ρ::Real
    Δt::Real
    pidx::UnitRange{Int}
    vidx::UnitRange{Int}
    guard::Function
    reset::Function
    decoupled_Φ::Function
    nominal::SafetyFilter
end

"""
"""
function(filt::ImpactAwareSafetyFilter)(
    Φargs::Tuple,
    Φ̇args::Tuple,
    x::Vector{<:Real},
    uref::Vector{<:Real};
    verbose::Bool = false
)::Vector{<:Real}
    if x[filt.vidx][end] < 0
        f = reshape(filt.nominal.flow.unactuated(x), length(x))
        g = reshape(filt.nominal.flow.actuated(x), length(x), filt.nominal.nu)
        x1 = x + filt.Δt*(f + g*filt.nominal.uvar)

        #vJ = min(1.0, exp(-filt.q*filt.guard(x))) * filt.reset(x1)[filt.vidx]
        vJ = filt.reset(x1)[filt.vidx]
        pJ = x[filt.pidx]
        Φ = filt.decoupled_Φ(Φargs..., pJ, vJ)
        #@show Φ
        @constraint(filt.nominal.model, Φ .+ filt.ρ .<= 0.0)
    end
    return filt.nominal(Φargs, Φ̇args, x, uref; verbose)
end
