/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.Actions
import SPG.Core.Algebra.Group
import SPG.Core.Geometry.SpatialOps
import SPG.Core.Geometry.SpinOps
import Mathlib.Data.Matrix.Basic

namespace SPG.Model.Multiferroic

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps

/--
Check if a group allows a non-zero Ferromagnetic moment (Axial Vector M).
Returns true if there exists a direction v such that for all g in group, g v = v.
This is equivalent to finding the invariant subspace of the axial vector representation.
For simplicity, we check if any basis vector (x, y, z) or their sum is invariant.
-/
def allows_ferromagnetism (group : List SPGElement) : List String :=
  let basis := [
    (![1, 0, 0], "Mx"),
    (![0, 1, 0], "My"),
    (![0, 0, 1], "Mz")
  ]
  basis.filterMap fun (v, name) =>
    if group.all fun g => act_on_axial_vector g v == v then
      some name
    else
      none

/--
Check if a group allows a non-zero Ferroelectric polarization (Polar Vector P).
-/
def allows_ferroelectricity (group : List SPGElement) : List String :=
  let basis := [
    (![1, 0, 0], "Px"),
    (![0, 1, 0], "Py"),
    (![0, 0, 1], "Pz")
  ]
  basis.filterMap fun (v, name) =>
    if group.all fun g => act_on_polar_vector g v == v then
      some name
    else
      none

/--
Check for allowed Magnetoelectric coupling terms (alpha_ij).
P_i = alpha_ij M_j
We check if alpha_ij is invariant under the group operations.
alpha' = g alpha
We need alpha'_ij == alpha_ij for all g.
The transformation rule for alpha (tensor connecting M (in) to P (out)) is:
act_on_tensor_rank2 g alpha .Axial .Polar
-/
def allowed_magnetoelectric_terms (group : List SPGElement) : List String :=
  let basis_indices := [(0, "x"), (1, "y"), (2, "z")]
  let terms := basis_indices.flatMap fun (i, i_name) =>
               basis_indices.map fun (j, j_name) =>
               ((i, j), s!"α_{i_name}{j_name}")

  -- We explicitly annotate the return type of filterMap to be List String
  terms.filterMap (fun (arg : ((ℕ × ℕ) × String)) =>
    match arg with
    | ((i, j), name) =>
      let T_basis : Matrix (Fin 3) (Fin 3) ℚ := Matrix.of fun r c =>
        if r.val == i && c.val == j then 1 else 0

      let projected := group.foldl (fun acc g =>
          acc + act_on_tensor_rank2 g T_basis .Axial .Polar
        ) 0

      let i_fin : Fin 3 := match i with | 0 => 0 | 1 => 1 | 2 => 2 | _ => 0
      let j_fin : Fin 3 := match j with | 0 => 0 | 1 => 1 | 2 => 2 | _ => 0

      let val := projected i_fin j_fin
      if val != 0 then
        some name
      else
        (none : Option String)
  )

/--
Compute the invariant Magnetoelectric Tensor form.
Returns a list of allowed non-zero components and their relations.
For now, just returns list of allowed indices (i,j) where alpha_ij can be non-zero.
-/
def get_allowed_me_components (group : List SPGElement) : List String :=
  let basis_indices := [(0, "x"), (1, "y"), (2, "z")]

  let terms := basis_indices.flatMap fun (i, i_name) =>
               basis_indices.map fun (j, j_name) =>
               ((i, j), s!"α_{i_name}{j_name}")

  -- Explicitly help type inference by iterating and filtering manually if filterMap is tricky
  -- or just define the function separately.
  -- Let's use List.filterMap with explicit type annotation on the lambda.
  terms.filterMap (fun (arg : ((ℕ × ℕ) × String)) =>
    match arg with
    | ((i, j), name) =>
      let T_basis : Matrix (Fin 3) (Fin 3) ℚ := Matrix.of fun r c =>
        if r.val == i && c.val == j then 1 else 0

      let projected := group.foldl (fun acc g =>
          acc + act_on_tensor_rank2 g T_basis .Axial .Polar
        ) 0

      let i_fin : Fin 3 := match i with | 0 => 0 | 1 => 1 | 2 => 2 | _ => 0
      let j_fin : Fin 3 := match j with | 0 => 0 | 1 => 1 | 2 => 2 | _ => 0

      let val := projected i_fin j_fin
      if val != 0 then
        some name
      else
        (none : Option String)
  )

end SPG.Model.Multiferroic
