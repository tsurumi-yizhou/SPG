/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Algebra.Actions
import Mathlib.Data.Matrix.Basic

namespace SPG.Model.Coupling

open SPG.Core.Algebra

/--
A generic Order Parameter type description.
Currently supports Polar Vector (P, u, E) and Axial Vector (M, S, H, L).
-/
inductive OrderType
| Polar
| Axial
deriving DecidableEq, Repr

/--
Check if a Trilinear Coupling term is invariant under the group.
Term form: C * O1_i * O2_j * O3_k
where O1, O2, O3 are vector order parameters.
Indices i, j, k are components (0=x, 1=y, 2=z).
-/
def check_trilinear_invariant (group : List SPGElement)
                              (t1 t2 t3 : OrderType)
                              (i j k : Fin 3) : Bool :=
  -- To check invariance of the term T = v1_i * v2_j * v3_k:
  -- We form the tensor product basis element E_ijk.
  -- We transform it: g(E_ijk) = sum_{pqr} R1_ip R2_jq R3_kr E_pqr
  -- And we check if the projection onto the identity representation is non-zero at component (i,j,k).
  -- Projector P = sum_g g.
  -- If (sum_g g(E_ijk)) has a non-zero component at (i,j,k), then this term is allowed (or coupled to itself).
  -- Actually, we want to know if the scalar constructed from these fields is invariant.
  -- The scalar is S = C * v1_i * v2_j * v3_k.
  -- g(S) = C * (g v1)_i * (g v2)_j * (g v3)_k ?
  -- No. If we write the energy as E = alpha * v1_i * v2_j * v3_k,
  -- We need E to be invariant.
  -- This means alpha_ijk must be an invariant tensor.
  -- So we check if the (i,j,k) component of the invariant tensor projector is non-zero.

  let projected_val := group.foldl (fun acc g =>
    -- Compute the (i,j,k) component of g(E_ijk)
    -- g(E_ijk) = (g v1) (times) (g v2) (times) (g v3)
    -- The (i,j,k) component of this tensor product is just:
    -- (g v1)_i * (g v2)_j * (g v3)_k ???
    -- NO.
    -- Let T = e_i (x) e_j (x) e_k.
    -- g T = (g e_i) (x) (g e_j) (x) (g e_k).
    -- We want to know the coefficient of e_i (x) e_j (x) e_k in the result sum_g g T.
    -- That coefficient is sum_g [ (g e_i)_i * (g e_j)_j * (g e_k)_k ].
    -- Let's verify.
    -- (g e_i) is the i-th column of the matrix representation of g acting on vector space 1.
    -- Let M1 be the matrix for g acting on v1. (g e_i) = sum_p M1_pi e_p.
    -- So (g e_i)_i is M1_ii.
    -- So we are summing M1_ii * M2_jj * M3_kk over g.
    -- This looks like checking if the diagonal element is non-zero.
    -- BUT, we need to be careful about basis.
    -- (g e_i) is a vector v1'. Its i-th component is v1'[i].
    -- So yes, we compute v1' = act(g, basis_i), v2' = act(g, basis_j), v3' = act(g, basis_k).
    -- Then we take product v1'[i] * v2'[j] * v3'[k].

    let mk_basis (idx : Fin 3) : Vec3 := fun n => if n == idx then 1 else 0

    let b1 := mk_basis i
    let b2 := mk_basis j
    let b3 := mk_basis k

    let b1_prime := match t1 with
      | .Polar => act_on_polar_vector g b1
      | .Axial => act_on_axial_vector g b1

    let b2_prime := match t2 with
      | .Polar => act_on_polar_vector g b2
      | .Axial => act_on_axial_vector g b2

    let b3_prime := match t3 with
      | .Polar => act_on_polar_vector g b3
      | .Axial => act_on_axial_vector g b3

    acc + (b1_prime i * b2_prime j * b3_prime k)
  ) (0 : ℚ)

  projected_val != 0

/--
Analyze Improper Ferroelectricity terms: P_i * M_j * M_k
(Coupling between Polarization and two Magnetization vectors)
Returns allowed non-zero components (i, j, k) where P_i M_j M_k is allowed.
-/
def analyze_improper_ferro (group : List SPGElement) : List String :=
  let indices := [0, 1, 2]
  let names := ["x", "y", "z"]

  indices.flatMap fun i =>
  indices.flatMap fun j =>
  indices.filterMap fun k =>
    let i_fin : Fin 3 := match i with | 0 => 0 | 1 => 1 | _ => 2
    let j_fin : Fin 3 := match j with | 0 => 0 | 1 => 1 | _ => 2
    let k_fin : Fin 3 := match k with | 0 => 0 | 1 => 1 | _ => 2

    if check_trilinear_invariant group .Polar .Axial .Axial i_fin j_fin k_fin then
      let name_i := names[i]!
      let name_j := names[j]!
      let name_k := names[k]!
      some s!"P{name_i} M{name_j} M{name_k}"
    else
      none

/--
Analyze DM-like interaction terms: P_i * L_j * (grad L)_k ?
Simplified: P_i * A_j * B_k where A, B are axial.
This is structurally identical to P M M if we consider L and grad L as just two axial vectors (or vector and tensor?).
Actually, the Lifshitz invariant is P . (L x curl L) - this involves derivatives.
In point group analysis, we usually look for P . (L x L') where L, L' are order parameters transforming like L.
This is exactly the same symmetry as P M M.
-/
def analyze_dm_like (group : List SPGElement) : List String :=
  analyze_improper_ferro group -- Same symmetry structure

end SPG.Model.Coupling
