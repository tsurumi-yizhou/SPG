/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Geometry.SpinOps
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.Module.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Algebra.Ring.Basic
import Mathlib.Algebra.Field.Defs
import Mathlib.Algebra.Ring.Rat
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace SPG.Core.Algebra

open SPG.Core.Geometry.SpinOps

/--
Action on a Polar Vector (e.g., Electric Polarization P, Displacement u).
Rules:
  - Spatial rotation R: v -> R v
  - Time reversal T: v -> v (Even)
-/
def act_on_polar_vector (g : SPGElement) (v : Vec3) : Vec3 :=
  Matrix.mulVec g.spatial v

/--
Action on an Axial Vector (e.g., Magnetization M, Spin S).
Rules:
  - Spatial rotation R: v -> (det R) * R v
  - Time reversal T: v -> -v (Odd)
-/
def act_on_axial_vector (g : SPGElement) (v : Vec3) : Vec3 :=
  let detR := Matrix.det g.spatial
  let spatial_transformed := detR • (Matrix.mulVec g.spatial v)
  if g.spin == spin_neg_I then -- Assuming spin_neg_I represents Time Reversal
    -spatial_transformed
  else
    spatial_transformed

/--
Action on a rank-2 Tensor connecting two vectors: Y_i = T_ij X_j
Transformation: T' = R T R^T (simplified, depends on X and Y types)

Case 1: Magnetoelectric Tensor alpha (P_i = alpha_ij M_j or P = alpha H)
P is Polar (Odd under I, Even under T)
H is Axial (Even under I, Odd under T)
alpha transforms to connect them.

General Tensor Action helper.
Transforms T_ij under g.
X_type: transformation of input vector X
Y_type: transformation of output vector Y
rule: Y' = g_Y Y, X' = g_X X
Y = T X => Y' = T' X'
g_Y Y = T' (g_X X) => Y = (g_Y^-1 T' g_X) X
So T = g_Y^-1 T' g_X  => T' = g_Y T g_X^-1

For orthogonal spatial rotations, g^-1 = g^T.
We need to handle the scalar factors (parity, time reversal) separately.
-/

inductive VectorType
| Polar
| Axial
deriving DecidableEq

def get_factor (g : SPGElement) (type : VectorType) : ℚ :=
  match type with
  | .Polar => (1 : ℚ)
  | .Axial =>
    let detR := Matrix.det g.spatial
    let time_rev := if g.spin == spin_neg_I then -1 else 1
    detR * time_rev

/--
Transform a matrix M (representing a tensor T_ij) under group element g.
input_type: The type of vector v_j in T_ij v_j
output_type: The type of vector w_i in w_i = T_ij v_j

Formula: T' = (factor_out * R) * T * (factor_in * R)^-1
            = (factor_out/factor_in) * R * T * R^T
-/
def act_on_tensor_rank2 (g : SPGElement) (T : Matrix (Fin 3) (Fin 3) ℚ) (in_type out_type : VectorType) : Matrix (Fin 3) (Fin 3) ℚ :=
  let f_in := get_factor g in_type
  let f_out := get_factor g out_type
  let factor := f_out / f_in -- division in Rationals

  -- R * T * R^T
  let R := g.spatial
  let RT := R.transpose
  let RTRT := R * T * RT

  factor • RTRT

end SPG.Core.Algebra
