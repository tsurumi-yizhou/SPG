/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Algebra.Group
import SPG.Core.Algebra.Actions
import SPG.Core.Geometry.SpatialOps
import SPG.Core.Geometry.SpinOps

namespace SPG.Model.Perturbation

open SPG.Core.Algebra

/--
Apply an Electric Field perturbation along a specific direction.
Reduces the symmetry group to the subgroup that leaves the E-field vector invariant.
Electric Field is a Polar Vector.
-/
def apply_electric_field (group : List SPGElement) (direction : Vec3) : List SPGElement :=
  group.filter fun g =>
    let E_prime := act_on_polar_vector g direction
    E_prime == direction

/--
Apply a Magnetic Field perturbation along a specific direction.
Magnetic Field is an Axial Vector.
-/
def apply_magnetic_field (group : List SPGElement) (direction : Vec3) : List SPGElement :=
  group.filter fun g =>
    let H_prime := act_on_axial_vector g direction
    H_prime == direction

/--
Apply Strain perturbation (simplified).
Strain is a symmetric rank-2 polar tensor.
For simplicity, we can model uniaxial strain along z as reducing symmetry like an E-field along z (for some groups)
or more precisely, preserving the tensor epsilon_zz.
Here we just provide a generic filter.
-/
def apply_strain (group : List SPGElement) (strain_tensor : Matrix (Fin 3) (Fin 3) ℚ) : List SPGElement :=
  group.filter fun g =>
    -- Strain transforms as epsilon' = R epsilon R^T (Polar Tensor)
    let epsilon_prime := act_on_tensor_rank2 g strain_tensor .Polar .Polar
    epsilon_prime == strain_tensor

end SPG.Model.Perturbation
