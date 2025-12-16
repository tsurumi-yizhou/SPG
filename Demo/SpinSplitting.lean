/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Algebra.Group
import SPG.Core.Geometry.SpatialOps
import SPG.Core.Geometry.SpinOps
import SPG.Interface.Notation
import SPG.Model.Hamiltonian
import SPG.Model.Perturbation
import Mathlib.Data.Matrix.Basic

namespace Demo.SpinSplitting

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps
open SPG.Interface
open SPG.Model.Hamiltonian
open SPG.Model.Perturbation

-- =================================================================
-- 1. Cubic Spin Splitting (PRL 2020)
-- =================================================================
-- Generators for Td (#216 F-43m)
def mat_4bar_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, 1, 0], ![-1, 0, 0], ![0, 0, -1]]
def mat_3_xyz  : Matrix (Fin 3) (Fin 3) ℚ := ![![0, 0, 1], ![1, 0, 0], ![0, 1, 0]]
def mat_m_xmy  : Matrix (Fin 3) (Fin 3) ℚ := ![![0, 1, 0], ![1, 0, 0], ![0, 0, 1]]

def gen_S4z : SPGElement := Op[mat_4bar_z, ^1]
def gen_C3  : SPGElement := Op[mat_3_xyz, ^1]
def gen_Mxmy : SPGElement := Op[mat_m_xmy, ^1]

def Td_Group : List SPGElement := generate_group [gen_S4z, gen_C3, gen_Mxmy]

def analyze_cubic : IO Unit := do
  IO.println "\n[1] Purely Cubic Spin Splittings (PRL 2020)"
  IO.println "    Symmetry: Td (Zinc Blende)"

  -- Check linear terms (should be zero for Td)
  let lin_terms := find_invariants Td_Group 1
  let has_linear := lin_terms.any (fun (p, s) => p.indices.length == 1 && s != .I)

  if has_linear then
    IO.println "    [UNEXPECTED] Found linear spin splitting!"
  else
    IO.println "    [OK] No linear spin splitting allowed (Dresselhaus linear forbidden)."

  -- Check cubic terms
  let cubic_invariants := find_invariants Td_Group 3
  let cubic_spin_terms := cubic_invariants.filter fun (p, s) => p.indices.length == 3 && s != .I

  IO.println s!"    Found {cubic_spin_terms.length} allowed cubic spin terms."
  IO.println "    Example allowed terms (indices):"
  for (p, s) in cubic_spin_terms.take 5 do
    IO.println s!"      {repr p.indices} * {repr s}"

  if cubic_spin_terms.length > 0 then
    IO.println "    => Cubic Dresselhaus effect is ALLOWED."

-- =================================================================
-- 2. E-field Controlled Zeeman Effect (PRL 2022)
-- =================================================================
-- Symmetry: Centrosymmetric Antiferromagnet.
-- Example: MnF2 (D4h).
-- D4h has inversion symmetry I.
-- Inversion I: k -> -k, sigma -> sigma.
-- Term k * sigma: (-k) * sigma = -(k * sigma). Forbidden!
-- So linear Zeeman effect is forbidden in D4h.
-- BUT, if we apply E-field, I is broken.

def mat_4_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, -1, 0], ![1, 0, 0], ![0, 0, 1]]
def mat_inv : Matrix (Fin 3) (Fin 3) ℚ := ![![ -1, 0, 0], ![0, -1, 0], ![0, 0, -1]]
def mat_2_x : Matrix (Fin 3) (Fin 3) ℚ := ![![1, 0, 0], ![0, -1, 0], ![0, 0, -1]]

-- Antiferromagnet with PT symmetry
-- Generators: C4z, PT, C2x (swaps sublattices)
def gen_AFM_C4z : SPGElement := Op[mat_4_z, ^1]
def gen_AFM_PT  : SPGElement := Op[mat_inv, ^-1] -- P * T
def gen_AFM_C2x : SPGElement := Op[mat_2_x, ^1]

def AFM_Group_D4h : List SPGElement := generate_group [gen_AFM_C4z, gen_AFM_PT, gen_AFM_C2x]

def analyze_efield_zeeman : IO Unit := do
  IO.println "\n[2] E-field Controlled Zeeman Effect (PRL 2022)"
  IO.println "    Material: Centrosymmetric AFM (D4h, PT-symmetric)"

  let group := AFM_Group_D4h

  -- 1. Check original symmetry
  let lin_orig := find_invariants group 1 |>.filter (fun (p, s) => p.indices.length == 1 && s != .I)
  if lin_orig.isEmpty then
    IO.println "    [Original] No linear spin splitting (k*sigma forbidden by PT or I)."
  else
    IO.println s!"    [Original] Found linear terms: {lin_orig.length}"

  -- 2. Apply Electric Field along z (breaks I, preserves C4z)
  let E_z : Vec3 := ![0, 0, 1]
  let group_E := apply_electric_field group E_z

  IO.println s!"    [With E||z] Reduced Group Order: {group_E.length} (from {group.length})"

  let lin_E := find_invariants group_E 1 |>.filter (fun (p, s) => p.indices.length == 1 && s != .I)

  if lin_E.isEmpty then
    IO.println "    [With E||z] Still no linear splitting."
  else
    IO.println s!"    [With E||z] Found {lin_E.length} allowed linear terms!"
    for (p, s) in lin_E do
      IO.println s!"      {repr p.indices} * {repr s}"
    IO.println "    => E-field induced Zeeman splitting is ALLOWED."

end Demo.SpinSplitting

def main : IO Unit := do
  Demo.SpinSplitting.analyze_cubic
  Demo.SpinSplitting.analyze_efield_zeeman
