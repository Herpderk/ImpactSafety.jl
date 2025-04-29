module ImpactSafety

using LinearAlgebra
using SparseArrays
using ForwardDiff
using JuMP
using Clarabel
using OSQP
using Plots

using HybridRobotDynamics:
        ControlAffineFlow,
        ManipulatorEquation,
        manipulator_inverses

export
        SafetyFilter,
        ImpactAwareSafetyFilter

include("impact.jl")
include("filter.jl")

end # module ImpactSafety
