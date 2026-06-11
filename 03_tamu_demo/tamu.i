# global parameters
rho = 1
mu = 3.333333e-3 # 1/Re = 1/3000
mesh = 'tamu.msh'
velocity_interp_method = 'rc'
advected_interp_method = 'upwind'

[GlobalParams]
  rhie_chow_user_object = ${velocity_interp_method}
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
  [mesh_file]
    type = FileMeshGenerator
    file = ${mesh}
  []
  [scale]
    type = TransformGenerator
    input = mesh_file
    transform = SCALE
    vector_value ='0.05249344 0.05249344 0.05249344'  # 1/d_in = 1/19.05
  []
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[Functions]
  [uz_in]
    type = ParsedFunction
    expression = -1*(60/49)*(1-y/0.5)^(1/7) # -u_avg*(60/49)*(1-r/R)^(1/7), with u_avg=1, R=0.5
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
    momentum_component = 'x'
    variable = 'uz'
    rho = ${rho}
  []
  [uz_advection]
    type = INSFVMomentumAdvection
    momentum_component = 'x'
    variable = uz
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
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
  [ur_time]
    type = INSFVMomentumTimeDerivative
    momentum_component = 'y'
    variable = ur
    rho = ${rho}
  []
  [ur_advection]
    type = INSFVMomentumAdvection
    momentum_component = 'y'
    variable = ur
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
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
  [ur_diffusion_rz]
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
[]

[AuxKernels]
  [mixing_len_aux_ker]
    type = WallDistanceMixingLengthAux
    walls = 'wall'
    variable = lm
    von_karman_const = 0.41
    delta = 0.5
    execute_on = 'initial'
  []
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
  [inlet-uz]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = uz
    functor = 'uz_in'
  []
  [inlet-ur]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = ur
    functor = 0
  []
  [no-slip-wall-uz]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = uz
    function = 0
  []
  [no-slip-wall-ur]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = ur
    function = 0
  []
  [outlet-p]
    type = INSFVOutletPressureBC
    boundary = 'outlet'
    variable = p
    function = 0
  []
  [axis-uz]
    type = INSFVSymmetryVelocityBC
    boundary = 'SYM'
    variable = uz
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = x
  []
  [axis-ur]
    type = INSFVSymmetryVelocityBC
    boundary = 'SYM'
    variable = ur
    u = uz
    v = ur
    mu = mu_eff
    momentum_component = y
  []
  [axis-p]
    type = INSFVSymmetryPressureBC
    boundary = 'SYM'
    variable = p
  []
[]

[VectorPostprocessors]
  [vpp]
    type = ElementValueSampler
    variable = 'uz ur p lm mu_t_out mu_eff_out'
    sort_by = x
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -ksp_gmres_restart'
    petsc_options_value = 'lu 100'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = IterationAdaptiveDT
    growth_factor = 1.25
    optimal_iterations = 8
    linear_iteration_ratio = 150
    dt = 0.5
    cutback_factor = 0.75
    cutback_factor_at_failure = 0.75
  []
  dtmin = 1e-6
  dtmax = 200
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 50
  l_tol = 1e-5
  l_max_its = 100
  start_time = 0
  end_time  = 10000
  num_steps = 10000
  steady_state_detection = true
  steady_state_tolerance = 1.e-6
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    execute_on = FINAL
    file_base = tamu
  []
  [csv]
    type = CSV
    execute_on = FINAL
    file_base = tamu
  []
  #[out]
  #  type = Checkpoint
  #  execute_on = FINAL
  #  file_base = tamu
  #[]
[]