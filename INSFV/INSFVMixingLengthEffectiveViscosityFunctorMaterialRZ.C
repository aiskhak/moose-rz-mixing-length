#include "INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.h"

registerMooseObject("NavierStokesApp", INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ);

InputParameters
INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ::validParams()
{
  InputParameters params = FunctorMaterial::validParams();

  params.addClassDescription(
    "RZ-aware effective dynamic viscosity functor material for the INSFV mixing-length "
    "model. The turbulent contribution is mu_t = rho * l_m^2 * |S|_RZ, where |S|_RZ "
    "is computed from the cylindrical-coordinate strain-rate invariant including the "
    "hoop-strain contribution u_r / r. The effective viscosity is mu_eff = mu + mu_t.");

  params.addRequiredParam<MooseFunctorName>(
      "property_name", "Name of the effective dynamic viscosity property.");

  params.addParam<MooseFunctorName>(
      "turbulent_viscosity_property_name",
      "Optional name of the turbulent dynamic viscosity property.");

  params.addRequiredParam<MooseFunctorName>("molecular_viscosity", "Molecular dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("rho", "Density.");
  params.addRequiredParam<MooseFunctorName>("mixing_length", "Mixing length.");

  params.addRequiredParam<MooseFunctorName>("u", "Velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "Velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "Velocity in the z direction.");

  params.set<ExecFlagEnum>("execute_on") = EXEC_ALWAYS;

  return params;
}

INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ::
    INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ(const InputParameters & params)
  : FunctorMaterial(params),
    _dim(blocksMaxDimension()),
    _rz_radial_coord(_subproblem.getAxisymmetricRadialCoord()),
    _molecular_viscosity(getFunctor<ADReal>("molecular_viscosity")),
    _rho(getFunctor<ADReal>("rho")),
    _mixing_length(getFunctor<ADReal>("mixing_length")),
    _u(getFunctor<ADReal>("u")),
    _v(params.isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _w(params.isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr)
{
  if (getBlockCoordSystem() != Moose::COORD_RZ)
    mooseError(name(), " is only valid on RZ blocks.");

  if (_dim >= 2 && !_v)
    mooseError("In two or more dimensions, the v velocity must be supplied using the 'v' parameter.");

  if (_dim >= 3 && !_w)
    mooseError("In three dimensions, the w velocity must be supplied using the 'w' parameter.");

  if (_rz_radial_coord >= _dim)
    mooseError(name(), ": the RZ radial coordinate is inconsistent with the block dimension.");

  const auto effective_property_name = getParam<MooseFunctorName>("property_name");

  addFunctorProperty<ADReal>(
      effective_property_name,
      [this](const auto & r, const auto & t) -> ADReal
      {
        return _molecular_viscosity(r, t) + computeTurbulentViscosity(r, t);
      });

  if (isParamValid("turbulent_viscosity_property_name"))
  {
    const auto turbulent_property_name =
        getParam<MooseFunctorName>("turbulent_viscosity_property_name");

    addFunctorProperty<ADReal>(
        turbulent_property_name,
        [this](const auto & r, const auto & t) -> ADReal
        {
          return computeTurbulentViscosity(r, t);
        });
  }
}