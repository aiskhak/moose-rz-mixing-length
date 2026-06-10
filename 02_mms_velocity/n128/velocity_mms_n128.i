
# 02_mms_velocity / n128
# Full RZ mixing-length solver MMS using final effective-viscosity material.

mu  = 1.0000000000000000e-02
rho = 1.0000000000000000e+00
lm0 = 2.0000000000000001e-01

advected_interp_method = 'average'
velocity_interp_method = 'rc'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = 128
    ny = 128
  []
  coord_type = 'RZ'
  rz_coord_axis = x # x=z; y=r
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[GlobalParams]
  rhie_chow_user_object = rc
[]

[Variables]
  [uz]
    type = INSFVVelocityVariable
    initial_condition = 0
  []

  [ur]
    type = INSFVVelocityVariable
    initial_condition = 0
  []

  [p]
    type = INSFVPressureVariable
    initial_condition = 0
  []

  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[AuxVariables]
  [lm]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []

  [mu_t_out]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []

  [mu_eff_out]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[ICs]
  [lm_ic]
    type = ConstantIC
    variable = lm
    value = ${lm0}
  []
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = uz
    v = ur
    pressure = p
  []
[]

[FunctorMaterials]
  [mixing_length_viscosity]
    type = INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ
    property_name = mu_eff
    turbulent_viscosity_property_name = mu_t
    molecular_viscosity = ${mu}
    rho = ${rho}
    mixing_length = lm
    u = uz
    v = ur
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = p
    rho = ${rho}
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
  []

  [pressure_pin]
    type = FVPointValueConstraint
    lambda = lambda
    variable = p
    point = '0.5 0.5 0'
    phi0 = 1.2500000000000000e-01
  []

  [uz_advection]
    type = INSFVMomentumAdvection
    variable = uz
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []

  [uz_diffusion]
    type = INSFVMomentumDiffusion
    variable = uz
    mu = mu_eff
    momentum_component = 'x'
    complete_expansion = true
    u = uz
    v = ur
  []

  [uz_pressure]
    type = INSFVMomentumPressure
    variable = uz
    momentum_component = 'x'
    pressure = p
  []

  [uz_body]
    type = INSFVBodyForce
    variable = uz
    functor = forcing_uz
    momentum_component = 'x'
  []

  [ur_advection]
    type = INSFVMomentumAdvection
    variable = ur
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []

  [ur_diffusion]
    type = INSFVMomentumDiffusion
    variable = ur
    mu = mu_eff
    momentum_component = 'y'
    complete_expansion = true
    u = uz
    v = ur
  []

  [ur_viscous_source_rz]
    type = INSFVMomentumViscousSourceRZ
    variable = ur
    mu = mu_eff
    momentum_component = 'y'
    complete_expansion = true
  []

  [ur_pressure]
    type = INSFVMomentumPressure
    variable = ur
    momentum_component = 'y'
    pressure = p
  []

  [ur_body]
    type = INSFVBodyForce
    variable = ur
    functor = forcing_ur
    momentum_component = 'y'
  []
[]

[AuxKernels]
  [mu_t_out_aux]
    type = FunctorAux
    variable = mu_t_out
    functor = mu_t
    execute_on = 'initial timestep_end final'
  []

  [mu_eff_out_aux]
    type = FunctorAux
    variable = mu_eff_out
    functor = mu_eff
    execute_on = 'initial timestep_end final'
  []
[]

[FVBCs]
  [axis_uz]
    type = INSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = uz
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = x
  []

  [axis_ur]
    type = INSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = ur
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = y
  []

  [axis_p]
    type = INSFVSymmetryPressureBC
    boundary = 'bottom'
    variable = p
  []

  [dirichlet_uz]
    type = INSFVNoSlipWallBC
    boundary = 'left right top'
    variable = uz
    function = exact_uz
  []

  [dirichlet_ur]
    type = INSFVNoSlipWallBC
    boundary = 'left right top'
    variable = ur
    function = exact_ur
  []
[]

[Functions]
  [exact_uz]
    type = ParsedFunction
    expression = 'x^2*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)'
  []

  [exact_ur]
    type = ParsedFunction
    expression = '(-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(y^4 - 2*y^3 + y^2)'
  []

  [exact_p]
    type = ParsedFunction
    expression = 'x*y^2'
  []

  [forcing_uz]
    type = ParsedFunction
    expression = '1.0*x^2*(1 - x)^2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(15*y^2 - 16*y + 3)*(y^4 - 2*y^3 + y^2) + 1.0*x^2*(1 - x)^2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))*(5*y^3 - 8*y^2 + 3*y) + 1.0*y^2 - 1.0*(x^2*(1 - x)^2*(30*y - 16) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y))*(0.04*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.01) - 0.04*(x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))*((-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(24*y^2 - 24*y + 4)*(4*y^3 - 6*y^2 + 2*y) + (2*x^2*(1 - x)^2*(30*y - 16) + 2*(-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y))*(x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))/2 + (2*x^2*(2*x - 2)*(15*y^2 - 16*y + 3) + 4*x*(1 - x)^2*(15*y^2 - 16*y + 3))*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)) + (-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(8*y^3 - 12*y^2 + 4*y)*(y^4 - 2*y^3 + y^2)/y^2 - 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^3)/sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) - 0.08*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))*((-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(-4*x^2 - 8*x*(2*x - 2) - 4*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))*(2*x^2*(2*x - 2)*(15*y^2 - 16*y + 3) + 4*x*(1 - x)^2*(15*y^2 - 16*y + 3) + 2*(12 - 24*x)*(y^4 - 2*y^3 + y^2))/2 + (x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))*(4*x^2*(5*y^3 - 8*y^2 + 3*y) + 8*x*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 4*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)) + (-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(-4*x^2 - 8*x*(2*x - 2) - 4*(1 - x)^2)*(y^4 - 2*y^3 + y^2)^2/y^2)/sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) - 1.0*(0.08*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.02)*(2*x^2*(5*y^3 - 8*y^2 + 3*y) + 4*x*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)) - 1.0*(x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))*(0.04*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.01)/y'
  []

  [forcing_ur]
    type = ParsedFunction
    expression = '1.0*x^2*(1 - x)^2*(-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(5*y^3 - 8*y^2 + 3*y)*(y^4 - 2*y^3 + y^2) + 2.0*x*y + 1.0*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)*(y^4 - 2*y^3 + y^2) - 1.0*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(0.08*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.02)*(12*y^2 - 12*y + 2) - 0.08*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y)*((-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(24*y^2 - 24*y + 4)*(4*y^3 - 6*y^2 + 2*y) + (2*x^2*(1 - x)^2*(30*y - 16) + 2*(-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y))*(x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))/2 + (2*x^2*(2*x - 2)*(15*y^2 - 16*y + 3) + 4*x*(1 - x)^2*(15*y^2 - 16*y + 3))*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)) + (-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(8*y^3 - 12*y^2 + 4*y)*(y^4 - 2*y^3 + y^2)/y^2 - 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^3)/sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) - 0.04*(x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))*((-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(-4*x^2 - 8*x*(2*x - 2) - 4*(1 - x)^2)*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))*(2*x^2*(2*x - 2)*(15*y^2 - 16*y + 3) + 4*x*(1 - x)^2*(15*y^2 - 16*y + 3) + 2*(12 - 24*x)*(y^4 - 2*y^3 + y^2))/2 + (x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))*(4*x^2*(5*y^3 - 8*y^2 + 3*y) + 8*x*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 4*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y)) + (-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(-4*x^2 - 8*x*(2*x - 2) - 4*(1 - x)^2)*(y^4 - 2*y^3 + y^2)^2/y^2)/sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) - 1.0*(0.04*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.01)*(x^2*(2*x - 2)*(15*y^2 - 16*y + 3) + 2*x*(1 - x)^2*(15*y^2 - 16*y + 3) + (12 - 24*x)*(y^4 - 2*y^3 + y^2)) - 1.0*((-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(0.08*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.02)*(4*y^3 - 6*y^2 + 2*y) - (-x^2*(2*x - 2) - 2*x*(1 - x)^2)*(0.08*sqrt(2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(4*y^3 - 6*y^2 + 2*y)^2 + (x^2*(1 - x)^2*(15*y^2 - 16*y + 3) + (-2*x^2 - 4*x*(2*x - 2) - 2*(1 - x)^2)*(y^4 - 2*y^3 + y^2))^2 + 2*(x^2*(2*x - 2)*(5*y^3 - 8*y^2 + 3*y) + 2*x*(1 - x)^2*(5*y^3 - 8*y^2 + 3*y))^2 + 2*(-x^2*(2*x - 2) - 2*x*(1 - x)^2)^2*(y^4 - 2*y^3 + y^2)^2/y^2) + 0.02)*(y^4 - 2*y^3 + y^2)/y)/y'
  []
[]

[VectorPostprocessors]
  [vpp]
    type = ElementValueSampler
    variable = 'uz ur p mu_t_out mu_eff_out'
    sort_by = x
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'lu NONZERO 200'
  []
[]

[Executioner]
  type = Steady
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 50
  l_tol = 1e-8
  l_max_its = 200
[]

[Outputs]
  print_linear_residuals = false

  [exodus]
    type = Exodus
    execute_on = FINAL
    file_base = velocity_mms_n128_out
  []

  [csv]
    type = CSV
    execute_on = FINAL
    file_base = velocity_mms_n128_csv
  []
[]
