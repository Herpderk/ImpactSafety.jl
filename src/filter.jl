"""
"""
mutable struct SafetyFilter
    nu::Int
    ϵ::Real
    Φ::Function
    Φ̇ub::Function
    flow::ControlAffineFlow
    uvar::Vector{JuMP.VariableRef}
    model::JuMP.Model
    #settings::NamedTuple
end

function SafetyFilter(
    nu::Int,
    ϵ::Real,
    Φ::Function,
    Φ̇ub::Function,
    flow::ControlAffineFlow
)::SafetyFilter
    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.@variable(model, uvar[1:nu])
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
    JuMP.@objective(filt.model, Min, obj)

    # Set constraints
    for (i, ub) in enumerate(Φ̇ub)
        if !isinf(ub)
            JuMP.@constraint(filt.model, cexpr[i] <= ub + filt.ϵ)
        end
    end

    # Set solver options
    if verbose
        JuMP.set_optimizer_attribute(filt.model, "verbose", true)
    else
        JuMP.set_optimizer_attribute(filt.model, "verbose", false)
    end

    # Solve
    JuMP.optimize!(filt.model)
    u = JuMP.value.(filt.uvar)

    # Clear constraints
    for c in JuMP.all_constraints(filt.model; include_variable_in_set_constraints = false)
        JuMP.delete(filt.model, c)
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
        JuMP.@constraint(filt.nominal.model, Φ .+ filt.ρ .<= 0.0)
    end
    return filt.nominal(Φargs, Φ̇args, x, uref; verbose)
end


"""
"""
mutable struct mjSafetyFilter
    nu::Int
    ϵ::Real
    Φ::Function
    Φ̇ub::Function
    δΦδx::Function
    uvar::Vector{JuMP.VariableRef}
    optmodel::JuMP.Model
    #settings::NamedTuple
end

function mjSafetyFilter(
    nu::Int,
    ϵ::Real,
    Φ::Function,
    Φ̇ub::Function,
    δΦδx::Function,
)::mjSafetyFilter
    optmodel = JuMP.Model(OSQP.Optimizer)
    JuMP.@variable(optmodel, uvar[1:nu])
    return mjSafetyFilter(nu, ϵ, Φ, Φ̇ub, δΦδx, uvar, optmodel)
end

"""
"""
function (filt::mjSafetyFilter)(
    Φargs::Tuple,
    Φ̇args::Tuple,
    model::MuJoCo.Model,
    data::MuJoCo.Data,
    bodyidx::Int;
    verbose::Bool = false
)::Nothing
    # Compute safe control set bounds
    Φ̇ub = filt.Φ̇ub(Φ̇args..., model, data, bodyidx)
    δΦδx = filt.δΦδx(Φargs, model, data, bodyidx)

    # Compute affine dynamics terms
    MuJoCo.mj_forward(model, data)
    f = [data.qvel; data.qacc - data.qfrc_actuator]
    g = [zeros(model.nv, model.nu); data.actuator_moment']
    cexpr = δΦδx * (f + g*filt.uvar)

    # Set objective
    uref = reshape(data.ctrl, length(data.ctrl))
    obj = filt.uvar'*filt.uvar - 2*uref'*filt.uvar + uref'*uref
    JuMP.@objective(filt.optmodel, Min, obj)

    # Set constraints
    for (i, ub) in enumerate(Φ̇ub)
        if !isinf(ub)
            JuMP.@constraint(filt.optmodel, cexpr[i] <= ub + filt.ϵ)
        end
    end

    # Set solver options
    if verbose
        JuMP.set_optimizer_attribute(filt.optmodel, "verbose", true)
    else
        JuMP.set_optimizer_attribute(filt.optmodel, "verbose", false)
    end

    # Solve
    JuMP.optimize!(filt.optmodel)
    data.ctrl .= JuMP.value.(filt.uvar)

    # Clear constraints
    for c in JuMP.all_constraints(filt.optmodel; include_variable_in_set_constraints = false)
        JuMP.delete(filt.optmodel, c)
    end
    return nothing
end


"""
"""
mutable struct mjVelocityExplicitFilter
    ρ::Real
    decoupled_Φ::Function
    nominal::mjSafetyFilter
end

"""
"""
function(filt::mjVelocityExplicitFilter)(
    Φargs::Tuple,
    Φ̇args::Tuple,
    model::MuJoCo.Model,
    data::MuJoCo.Data,
    bodyidx::Int;
    ecoeff::Real = 0.1,
    verbose::Bool = false
)::Nothing
    if data.cvel[bodyidx, 3] < 0
        # Task state jacobian w.r.t. joint position
        Jtask = MuJoCo.mj_zeros(3, model.nv)
        MuJoCo.mj_jacBodyCom(model, data, Jtask, nothing, bodyidx)

        # Get velocity integration law
        fv = data.qacc - data.qfrc_actuator
        gv = data.actuator_moment'
        v1 = data.qvel + model.opt.timestep*(fv + gv*filt.nominal.uvar)

        # Get task space velocity and position
        vtask = -ecoeff *  Jtask * v1
        ptask = data.xipos[bodyidx, :]

        Φ = filt.decoupled_Φ(Φargs..., ptask, vtask)
        JuMP.@constraint(filt.nominal.optmodel, Φ .+ filt.ρ .<= 0.0)
    end
    return filt.nominal(Φargs, Φ̇args, model, data, bodyidx; verbose)
end
