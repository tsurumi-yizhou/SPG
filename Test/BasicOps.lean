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

open SPG
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface

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
  let op := Op[mat_2_x, ^-1]
  let v : Vec3 := ![1, 1, 1]
  let expected : Vec3 := ![-1, 1, 1]

  let res := magnetic_action op v
  if res == expected then
    IO.println "[PASS] Combined C2x * T action correct."
  else
    IO.println s!"[FAIL] Combined action wrong. Expected {repr expected}, got {repr res}."

def main : IO Unit := do
  IO.println "=== Running Basic Operations Tests ==="
  test_identity_action
  test_time_reversal
  test_combined_op
  IO.println "======================================"
