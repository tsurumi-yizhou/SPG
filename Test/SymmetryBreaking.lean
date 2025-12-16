/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import SPG.Algebra.Group
import SPG.Physics.SymmetryBreaking
import SPG.Geometry.SpatialOps
import SPG.Geometry.SpinOps
import SPG.Interface.Notation

open SPG
open SPG.Physics
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface
open SPG.Algebra

-- Test Scenario: Ferromagnetic Iron (bcc Fe) model
-- Point group Oh (m3m).
-- Let's take a subgroup D4h (4/mmm) for simplicity (tetragonal).
-- Generators: C4z, C2x, I.
-- If we apply FM order along z, what is the MPG?
-- FM along z means z-vector is invariant.
-- C4z keeps z invariant. C2x flips z -> -z. I flips z -> -z.
-- So C4z is in MPG. C2x is NOT. I is NOT.
-- Wait, magnetic moment is an axial vector (pseudovector).
-- Under Inversion: r -> -r, but magnetic moment m -> m (axial).
-- So I keeps m invariant!
-- Under C2x (rotation 180 about x): y->-y, z->-z. Axial vector components:
-- mx -> mx, my -> -my, mz -> -mz.
-- So C2x flips mz. Not in MPG.
-- So MPG should contain C4z and I. (C4h point group).

def mat_4_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, -1, 0], ![1, 0, 0], ![0, 0, 1]]
def mat_inv : Matrix (Fin 3) (Fin 3) ℚ := ![![ -1, 0, 0], ![0, -1, 0], ![0, 0, -1]]

def gen_C4z : SPGElement := Op[mat_4_z, ^1]
def gen_C2x : SPGElement := Op[mat_2_x, ^1]
def gen_Inv : SPGElement := Op[mat_inv, ^1]

def full_group : List SPGElement := generate_group [gen_C4z, gen_C2x, gen_Inv]

def test_ferromagnet_z : IO Unit := do
  let neel_z : Vec3 := ![0, 0, 1]
  let mpg := get_mpg full_group neel_z

  -- We expect C4z and Inv to be in MPG.
  -- C2x should NOT be in MPG.

  let has_c4z := mpg.contains gen_C4z
  let has_inv := mpg.contains gen_Inv
  let has_c2x := mpg.contains gen_C2x

  if has_c4z && has_inv && !has_c2x then
    IO.println "[PASS] Ferromagnet (z-axis) MPG contains C4z and Inv, excludes C2x."
  else
    IO.println s!"[FAIL] FM MPG incorrect. C4z:{has_c4z}, Inv:{has_inv}, C2x:{has_c2x}."

-- Test Scenario: Antiferromagnet / Altermagnet logic
-- If we have C2x * T (Time Reversal).
-- C2x flips mz -> -mz. T flips mz -> -mz.
-- So C2x * T should keep mz invariant!
def test_afm_symmetry : IO Unit := do
  let op_afm : SPGElement := Op[mat_2_x, ^-1]
  let neel_z : Vec3 := ![0, 0, 1]

  -- Check if op_afm preserves neel_z
  let v_prime := magnetic_action op_afm neel_z

  if v_prime == neel_z then
    IO.println "[PASS] C2x * T preserves z-moment (AFM symmetry)."
  else
    IO.println s!"[FAIL] C2x * T should preserve z-moment. Got {repr v_prime}."

def main : IO Unit := do
  IO.println "=== Running Symmetry Breaking Tests ==="
  test_ferromagnet_z
  test_afm_symmetry
  IO.println "======================================="
