module ImpactSafety

using LinearAlgebra
using SparseArrays
using ForwardDiff
import JuMP
using Clarabel
using OSQP
using Plots
import MuJoCo

using HybridRobotDynamics:
        ControlAffineFlow,
        ManipulatorEquation,
        manipulator_inverses

export
        SafetyFilter,
        ImpactAwareSafetyFilter,
        mjSafetyFilter,
        mjVelocityExplicitFilter,
        taskspace_Jdot

include("impact.jl")
include("filter.jl")
include("finitediff.jl")

end # module ImpactSafety
