#pragma once

#include "FunctorMaterial.h"
#include "MooseTypes.h"
#include "metaphysicl/raw_type.h"

class INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ : public FunctorMaterial
{
public:
  static InputParameters validParams();

  INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ(const InputParameters & parameters);

protected:
  template <typename SpaceArg, typename TimeArg>
  ADReal computeTurbulentViscosity(const SpaceArg & r, const TimeArg & t) const;

  const unsigned int _dim;
  const unsigned int _rz_radial_coord;

  const Moose::Functor<ADReal> & _molecular_viscosity;
  const Moose::Functor<ADReal> & _rho;
  const Moose::Functor<ADReal> & _mixing_length;

  const Moose::Functor<ADReal> & _u;
  const Moose::Functor<ADReal> * const _v;
  const Moose::Functor<ADReal> * const _w;
};

template <typename SpaceArg, typename TimeArg>
ADReal
INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ::computeTurbulentViscosity(
    const SpaceArg & r, const TimeArg & t) const
{
  using std::sqrt;

  //constexpr Real offset = 1e-15;
  constexpr Real radius_tol = 1e-14;

  const auto grad_u = _u.gradient(r, t);

  ADReal strain_norm_sq = 2.0 * Utility::pow<2>(grad_u(0));

  ADReal radial_velocity = 0.0;
  ADReal radial_velocity_radial_gradient = 0.0;

  if (_rz_radial_coord == 0)
  {
    radial_velocity = _u(r, t);
    radial_velocity_radial_gradient = grad_u(0);
  }

  if (_dim >= 2)
  {
    const auto grad_v = _v->gradient(r, t);

    strain_norm_sq +=
        2.0 * Utility::pow<2>(grad_v(1)) + Utility::pow<2>(grad_v(0) + grad_u(1));

    if (_rz_radial_coord == 1)
    {
      radial_velocity = (*_v)(r, t);
      radial_velocity_radial_gradient = grad_v(1);
    }

    if (_dim >= 3)
    {
      const auto grad_w = _w->gradient(r, t);

      strain_norm_sq +=
          2.0 * Utility::pow<2>(grad_w(2)) +
          Utility::pow<2>(grad_u(2) + grad_w(0)) +
          Utility::pow<2>(grad_v(2) + grad_w(1));

      if (_rz_radial_coord == 2)
      {
        radial_velocity = (*_w)(r, t);
        radial_velocity_radial_gradient = grad_w(2);
      }
    }
  }

  ADReal radius = 1.0;

  if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::ElemArg>)
    radius = r.elem->vertex_average()(_rz_radial_coord);
  else if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::FaceArg>)
    radius = r.fi->faceCentroid()(_rz_radial_coord);
  else if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::ElemQpArg>)
    radius = r.elem->vertex_average()(_rz_radial_coord);
  else if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::ElemSideQpArg>)
    radius = r.elem->vertex_average()(_rz_radial_coord);
  else if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::ElemPointArg>)
    radius = r.elem->vertex_average()(_rz_radial_coord);
  else if constexpr (std::is_same_v<std::decay_t<SpaceArg>, Moose::NodeArg>)
    mooseError(name(),
               ": nodal evaluation is not supported for ",
               "INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ. ",
               "Use element, face, element-QP, element-side-QP, or element-point evaluation.");
  else
    mooseError(name(),
               ": unsupported functor evaluation type in ",
               "INSFVMixingLengthEffectiveViscosityFunctorMaterialRZ.");

  /*
   * In RZ, the hoop strain is u_r / r. At the axis, regularity gives
   * lim_{r -> 0} u_r / r = d u_r / d r. This avoids a 0/0 evaluation
   * when the functor is evaluated on an axis boundary face.
   */
  const ADReal hoop_strain =
      MetaPhysicL::raw_value(radius) > radius_tol ? radial_velocity / radius
                                                  : radial_velocity_radial_gradient;

  strain_norm_sq += 2.0 * Utility::pow<2>(hoop_strain);

  /*
   * The strain-rate magnitude is exactly zero when strain_norm_sq is zero.
   * We avoid evaluating the AD derivative of sqrt(x) at x = 0, where it is
   * singular. This is not a strain-rate regularization: the returned value is
   * still exactly zero for zero strain.
   */
  ADReal strain_norm = 0.0;
  if (MetaPhysicL::raw_value(strain_norm_sq) > 0.0)
    strain_norm = sqrt(strain_norm_sq);

  const ADReal lm = _mixing_length(r, t);

  return _rho(r, t) * lm * lm * strain_norm;
}