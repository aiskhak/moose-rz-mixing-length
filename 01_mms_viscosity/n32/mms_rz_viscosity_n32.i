
# 01_mms_viscosity / n32
#
# Grid-convergence verification of the RZ mixing-length eddy-viscosity
# evaluation using the final effective-viscosity functor material.
#
# This input does not solve a flow problem. It initializes uz and ur,
# evaluates mu_t and mu_eff from the solve-path functor material, and writes
# cell-centered values.
#
# Prescribed velocity:
#   uz = sin(pi*x)*(1+r^2)
#   ur = cos(pi*x)*r^2
#
# Domain:
#   x in [0.0, 1.0]
#   r in [0.5, 1.5]

rho = 1.0000000000000000e+00
mu  = 3.3333330000000000e-03
lm0 = 5.0000000000000000e-01

advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = uz
    v = ur
    pressure = p
  []
[]

[Mesh]
  coord_type = 'RZ'
  rz_coord_axis = x # x=z; y=r

  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 32
    ny = 32
    xmin = 0.0
    xmax = 1.0
    ymin = 0.5
    ymax = 1.5
  []
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[Variables]
  [uz]
    type = INSFVVelocityVariable
  []

  [ur]
    type = INSFVVelocityVariable
  []

  [p]
    type = INSFVPressureVariable
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

[Functions]
  [uz_exact]
    type = ParsedFunction
    expression = 'sin(pi*x)*(1+y^2)'
  []

  [ur_exact]
    type = ParsedFunction
    expression = 'cos(pi*x)*y^2'
  []
[]

[ICs]
  [uz_ic]
    type = FunctionIC
    variable = uz
    function = uz_exact
  []

  [ur_ic]
    type = FunctionIC
    variable = ur
    function = ur_exact
  []

  [p_ic]
    type = ConstantIC
    variable = p
    value = 0.0
  []

  [lm_ic]
    type = ConstantIC
    variable = lm
    value = ${lm0}
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
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []

  [uz_time]
    type = INSFVMomentumTimeDerivative
    variable = uz
    momentum_component = 'x'
    rho = ${rho}
  []

  [ur_time]
    type = INSFVMomentumTimeDerivative
    variable = ur
    momentum_component = 'y'
    rho = ${rho}
  []
[]


[AuxKernels]
  [mu_t_out_aux]
    type = FunctorAux
    variable = mu_t_out
    functor = mu_t
    execute_on = 'initial final'
  []

  [mu_eff_out_aux]
    type = FunctorAux
    variable = mu_eff_out
    functor = mu_eff
    execute_on = 'initial final'
  []
[]

[VectorPostprocessors]
  [vpp]
    type = ElementValueSampler
    variable = 'uz ur lm mu_t_out mu_eff_out'
    sort_by = x
    execute_on = 'initial final'
  []
[]

[Executioner]
  type = Transient
  num_steps = 0
[]

[Outputs]
  print_linear_residuals = false

  [csv]
    type = CSV
    execute_on = 'initial final'
    file_base = mms_rz_viscosity_n32_csv
  []

  [exodus]
    type = Exodus
    execute_on = 'initial final'
    file_base = mms_rz_viscosity_n32_out
  []
[]
