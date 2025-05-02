#= """
"""
function taskspace_Jdot(
    model::MuJoCo.Model,
    data::MuJoCo.Data,
    bodyidx::Int;
    epsilon::Real = 1e-12
)::Matrix
    # Save original config
    orig_q = copy(data.qpos)

    nq = model.nq
    nv = model.nv

    # Update forward dynamics
    MuJoCo.mj_forward(model, data)

    # Compute base Jacobian: J(q)
    jacp_0 = MuJoCo.mj_zeros(3, nv)
    MuJoCo.mj_jacBodyCom(model, data, jacp_0, nothing, bodyidx)
    xdot_0 = jacp_0 * data.qvel

    # Output: ∂ẋ/∂q ∈ ℝ^{3×nq}
    Jdot = MuJoCo.mj_zeros(3, nq)

    for k in 1:nq
        # Get perturbed config
        dq = MuJoCo.zeros(nq)
        dq[k] = epsilon
        data.qpos .+= dq
        MuJoCo.mj_forward(model, data)

        # Compute perturbed Jacobian
        jacp_k = MuJoCo.mj_zeros(3, nv)
        MuJoCo.mj_jacBodyCom(model, data, jacp_k, nothing, bodyidx)
        xdot_k = jacp_k * data.qvel
        Jdot[:, k] = (xdot_k - xdot_0) / epsilon
    end

    # Reset state
    data.qpos .= orig_q
    MuJoCo.mj_forward(model, data)
    return Jdot
end =#

function taskspace_Jdot(
    model::MuJoCo.Model,
    data::MuJoCo.Data,
    bodyidx::Int;
    δt::Real = 1e-6
)::Matrix
    # Backup original state
    qpos0 = copy(data.qpos)
    qvel0 = copy(data.qvel)
    Δt = model.opt.timestep

    # J(q)
    jacp0 = MuJoCo.mj_zeros(3, model.nv)
    MuJoCo.mj_forward(model, data)
    MuJoCo.mj_jacBodyCom(model, data, jacp0, nothing, bodyidx)

    # Step forward in time
    model.opt.timestep = δt
    MuJoCo.mj_step(model, data)

    # J(q + Δq)
    jacp1 = MuJoCo.mj_zeros(3, model.nv)
    MuJoCo.mj_forward(model, data)
    MuJoCo.mj_jacBodyCom(model, data, jacp1, nothing, bodyidx)

    # Restore state
    data.qpos .= qpos0
    data.qvel .= qvel0
    model.opt.timestep = Δt
    MuJoCo.mj_forward(model, data)

    # Compute finite difference: dJ/dt ≈ (J1 - J0) / Δt
    Jdot = (jacp1 - jacp0) / δt
    return Jdot
end

