/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import SPG.Algebra.Group
import SPG.Geometry.SpatialOps
import SPG.Geometry.SpinOps
import SPG.Interface.Notation
import SPG.Physics.Hamiltonian
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace Demo.Altermagnet

open SPG
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface
open SPG.Algebra
open SPG.Physics.Hamiltonian

-- 1. Generators for D4h altermagnet
def mat_4_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, -1, 0], ![1, 0, 0], ![0, 0, 1]]
def mat_inv : Matrix (Fin 3) (Fin 3) ℚ := ![![ -1, 0, 0], ![0, -1, 0], ![0, 0, -1]]
-- mat_2_xy defined explicitly for clarity (rotation by 180 degrees about x+y axis)
def mat_2_xy : Matrix (Fin 3) (Fin 3) ℚ := ![![0, 1, 0], ![1, 0, 0], ![0, 0, -1]]

-- Generators
-- C4z * T : The key altermagnetic symmetry
-- C2xy    : Rotation about x+y axis
-- I       : Inversion
def gen_C4z_TR : SPGElement := Op[mat_4_z, ^-1]
def gen_C2xy   : SPGElement := Op[mat_2_xy, ^1]
def gen_Inv    : SPGElement := Op[mat_inv, ^1]

def Altermagnet_Group : List SPGElement := generate_group [gen_C4z_TR, gen_C2xy, gen_Inv]

-- 2. ICE Symbol Output (Simplistic)
-- We just print generators for now, as full ICE symbol logic is complex.
def ice_symbol : String :=
  "4'22 (D4h magnetic)" -- Placeholder, deriving from generators

end Demo.Altermagnet

def main : IO Unit := do
  IO.println s!"Generated Group Size: {Demo.Altermagnet.Altermagnet_Group.length}"
  IO.println s!"ICE Symbol (Approx): {Demo.Altermagnet.ice_symbol}"
  IO.println "Allowed Hamiltonian Terms H(k) (up to quadratic):"

  let invariants := SPG.Physics.Hamiltonian.find_invariants Demo.Altermagnet.Altermagnet_Group
  for (p, s) in invariants do
    let p_str := match p with
      | .const => "1"
      | .x => "kx" | .y => "ky" | .z => "kz"
      | .xx => "kx^2" | .yy => "ky^2" | .zz => "kz^2"
      | .xy => "kx ky" | .yz => "ky kz" | .zx => "kz kx"
    let s_str := match s with
      | .I => "I" | .x => "σx" | .y => "σy" | .z => "σz"
    IO.println s!"  {p_str} * {s_str}"

  -- Manually check d-wave altermagnet term: (kx^2 - ky^2) * sigma_z?
  -- Or kx ky * sigma_z ?
  -- Let's define a custom checker for kx^2 - ky^2
  let check_custom (f : SPG.Vec3 → ℚ) (s : SPG.Physics.Hamiltonian.SpinComp) : Bool :=
    let test_ks : List SPG.Vec3 := [![1, 0, 0], ![0, 1, 0], ![0, 0, 1], ![1, 1, 0], ![1, 0, 1], ![0, 1, 1], ![1, 2, 3]]
    test_ks.all fun k =>
      -- Re-implement loop
      Demo.Altermagnet.Altermagnet_Group.all fun g =>
        let val_gk := f (SPG.Physics.Hamiltonian.act_on_k g k)
        let val_k  := f k
        -- Re-use projection logic from check_invariant
        if s == .I then
          val_gk == val_k
        else
          let s_prime := SPG.Physics.Hamiltonian.act_on_spin g s
          let coeff := SPG.Physics.Hamiltonian.project_spin s_prime s

          -- Check eigenstate property (no mixing)
          let is_eigen :=
              (s == .x && s_prime 1 == 0 && s_prime 2 == 0) ||
              (s == .y && s_prime 0 == 0 && s_prime 2 == 0) ||
              (s == .z && s_prime 0 == 0 && s_prime 1 == 0)

          if is_eigen then
             val_gk * coeff == val_k
          else
             false

  let kx2_minus_ky2 (k : SPG.Vec3) : ℚ := k 0 * k 0 - k 1 * k 1
  let kx2_plus_ky2  (k : SPG.Vec3) : ℚ := k 0 * k 0 + k 1 * k 1
  let kx_ky         (k : SPG.Vec3) : ℚ := k 0 * k 1

  if check_custom kx2_minus_ky2 .z then
    IO.println "  (kx^2 - ky^2) * σz  [d-wave altermagnetism (x^2-y^2 type)]"
  else
    IO.println "  (kx^2 - ky^2) * σz  [FORBIDDEN]"
    -- Analyze why it is forbidden (or allowed)
    let _ ← SPG.Physics.Hamiltonian.analyze_term_symmetry Demo.Altermagnet.Altermagnet_Group kx2_minus_ky2 .z "(kx^2 - ky^2)" "σz"

  if check_custom kx_ky .z then
    IO.println "  kx ky * σz          [d-wave altermagnetism (xy type)]"
  else
    IO.println "  kx ky * σz          [FORBIDDEN]"
    let _ ← SPG.Physics.Hamiltonian.analyze_term_symmetry Demo.Altermagnet.Altermagnet_Group kx_ky .z "kx ky" "σz"

  if check_custom kx2_plus_ky2 .I then
    IO.println "  (kx^2 + ky^2) * I   [Standard kinetic term]"
  else
    IO.println "  (kx^2 + ky^2) * I   [FORBIDDEN]"
    let _ ← SPG.Physics.Hamiltonian.analyze_term_symmetry Demo.Altermagnet.Altermagnet_Group kx2_plus_ky2 .I "(kx^2 + ky^2)" "I"
