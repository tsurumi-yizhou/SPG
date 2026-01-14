/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import SPG.Algebra.Actions
import SPG.Geometry.SpatialOps
import SPG.Geometry.SpinOps
import SPG.Interface.Notation
import SPG.Physics.Hamiltonian
import SPG.Data.MagneticGroups
import SPG.Data.Tetragonal

open SPG
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface
open SPG.Physics.Hamiltonian
open SPG.Data.MagneticGroups
open SPG.Data.Tetragonal

-- Test Case 1: Identity Action
-- Op[I, I] should leave any vector unchanged.
def test_identity_action : IO Unit := do
  let id_el := Op[1, ^1]
  let v : Vec3 := ![1, 2, 3]
  let v_prime := magnetic_action id_el v
  if v_prime == v then
    IO.println "[PASS] Identity action leaves vector unchanged."
  else
    IO.println s!"[FAIL] Identity action failed. Expected {repr v}, got {repr v_prime}."

-- Test Case 2: Time Reversal
-- Op[I, -I] should flip spin but not space (so magnetic action flips, electric action unchanged).
-- Note: magnetic_action applies BOTH spatial and spin matrices.
-- spin_neg_I is -1. So magnetic_action should return -v.
def test_time_reversal : IO Unit := do
  let tr_el := Op[1, ^-1]
  let v : Vec3 := ![1, 0, 0]

  let mag_v := magnetic_action tr_el v
  let elec_v := electric_action tr_el v

  if mag_v == -v then
    IO.println "[PASS] Time reversal flips magnetic vector."
  else
    IO.println s!"[FAIL] Time reversal magnetic action wrong. Expected {repr (-v)}, got {repr mag_v}."

  if elec_v == v then
    IO.println "[PASS] Time reversal leaves electric vector unchanged."
  else
    IO.println s!"[FAIL] Time reversal electric action wrong. Expected {repr v}, got {repr elec_v}."

-- Test Case 3: Combined Operation (C2x * T)
-- C2x: (x, y, z) -> (x, -y, -z)
-- T: flips sign
-- Combined: (x, y, z) -> (-x, y, z)
def test_combined_op : IO Unit := do
  let op := Op[SPG.Geometry.SpatialOps.mat_2_x, ^-1]
  let v : Vec3 := ![1, 1, 1]
  let expected : Vec3 := ![-1, 1, 1]

  let res := magnetic_action op v
  if res == expected then
    IO.println "[PASS] Combined C2x * T action correct."
  else
    IO.println s!"[FAIL] Combined action wrong. Expected {repr expected}, got {repr res}."

def expect_contains (xs : List String) (s : String) : Bool :=
  xs.any (fun t => t = s)

def test_invariant_solver_altermagnet : IO Unit := do
  let blocks := invariants_vector_by_degree_solve Altermagnet_Group_D4h 2
  let okAll :=
    blocks.all (fun (_, hs) => hs.all (fun h => isInvariantHam Altermagnet_Group_D4h h))
  let deg2 := blocks.find? (fun p => p.fst = 2)
  let okDeg2 :=
    match deg2 with
    | none => false
    | some (_, hs) =>
      let ss := hs.map ham_to_string
      expect_contains ss "(-ky^2 + kx^2)*σz"
  if okAll && okDeg2 then
    IO.println "[PASS] Invariant solver works for altermagnet."
  else
    IO.println "[FAIL] Invariant solver failed for altermagnet."

def test_invariant_solver_ferromagnet : IO Unit := do
  let blocks := invariants_vector_by_degree_solve Ferromagnet_Group_D4h_z 0
  let okAll :=
    blocks.all (fun (_, hs) => hs.all (fun h => isInvariantHam Ferromagnet_Group_D4h_z h))
  let deg0 := blocks.find? (fun p => p.fst = 0)
  let okDeg0 :=
    match deg0 with
    | none => false
    | some (_, hs) =>
      let ss := hs.map ham_to_string
      expect_contains ss "(1)*σz"
  if okAll && okDeg0 then
    IO.println "[PASS] Invariant solver works for ferromagnet."
  else
    IO.println "[FAIL] Invariant solver failed for ferromagnet."

def test_invariant_solver_altermagnet_complex : IO Unit := do
  let blocks := invariants_vector_by_degree_solveC Altermagnet_Group_D4h 1
  let okAll :=
    blocks.all (fun (_, hs) => hs.all (fun h => isInvariantCHam Altermagnet_Group_D4h h))
  if okAll then
    IO.println "[PASS] Complex invariant solver works for altermagnet."
  else
    IO.println "[FAIL] Complex invariant solver failed for altermagnet."

def test_invariant_solver_scalar_complex : IO Unit := do
  let blocks := invariants_scalar_by_degree_solveC Altermagnet_Group_D4h 1
  let okAll :=
    blocks.all (fun (_, ps) => ps.all (fun p => isInvariantCPoly Altermagnet_Group_D4h p))
  if okAll then
    IO.println "[PASS] Complex scalar invariant solver works for altermagnet."
  else
    IO.println "[FAIL] Complex scalar invariant solver failed for altermagnet."

def test_hermitian_projection : IO Unit := do
  let H : CKPHam := { scalar := ⟨1, 1⟩, vector := fun _ => 0 }
  let projected := project_hermitianCHam H
  let ok :=
    (!isHermitianCHam H) &&
      isHermitianCHam projected &&
      decide (projected.scalar = ⟨1, 0⟩)
  if ok then
    IO.println "[PASS] Hermitian projection produces Hermitian Hamiltonian."
  else
    IO.println "[FAIL] Hermitian projection test failed."

def test_combined_laue_magnetic_workflow : IO Unit := do
  let laue := Laue_D2d_gens
  let mag : List SPGElement := [SPG.Data.mk_ice_element mat_4bar_z true]
  let group := group_from_laue_magnetic laue mag
  let blocks := allowed_cham_by_degree_from_gens laue mag 1
  let okAll :=
    blocks.all (fun (_, hs) =>
      hs.all (fun h => isHermitianCHam h && isInvariantCHam group h)
    )
  if okAll then
    IO.println "[PASS] Combined Laue+magnetic workflow returns Hermitian invariants."
  else
    IO.println "[FAIL] Combined Laue+magnetic workflow test failed."

def main : IO Unit := do
  IO.println "=== Running Basic Operations Tests ==="
  test_identity_action
  test_time_reversal
  test_combined_op
  test_invariant_solver_altermagnet
  test_invariant_solver_ferromagnet
  test_invariant_solver_altermagnet_complex
  test_invariant_solver_scalar_complex
  test_hermitian_projection
  test_combined_laue_magnetic_workflow
  IO.println "======================================"
