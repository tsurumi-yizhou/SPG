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
import SPG.Material.MagneticGroup
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace Demo.MagneticPhases

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps
open SPG.Interface
open SPG.Model.Hamiltonian
open SPG.Material.MagneticGroup

def check_linear_k_sigma (group : List SPGElement) : IO Unit := do
  IO.println "  Checking linear spin-splitting terms (k * sigma):"
  let terms := [
    (Poly.x, SpinComp.x, "kx * σx"), (Poly.x, SpinComp.y, "kx * σy"), (Poly.x, SpinComp.z, "kx * σz"),
    (Poly.y, SpinComp.x, "ky * σx"), (Poly.y, SpinComp.y, "ky * σy"), (Poly.y, SpinComp.z, "ky * σz"),
    (Poly.z, SpinComp.x, "kz * σx"), (Poly.z, SpinComp.y, "kz * σy"), (Poly.z, SpinComp.z, "kz * σz")
  ]

  let mut found := false
  for (p, s, name) in terms do
    if check_invariant group p s then
      IO.println s!"    [ALLOWED] {name}"
      found := true

  if !found then
    IO.println "    [NONE] No linear spin-splitting terms found."

def analyze_phase (name : String) (group : List SPGElement) : IO Unit := do
  IO.println s!"\n========================================"
  IO.println s!"Phase Analysis: {name}"
  IO.println s!"Group Order: {group.length}"
  IO.println "----------------------------------------"

  -- Check Chemical Potential
  if check_invariant group Poly.const .I then
    IO.println "  [ALLOWED] Chemical Potential (1 * I)"
  else
    IO.println "  [FORBIDDEN] Chemical Potential (1 * I) - Strange!"

  -- Check Standard Kinetic Energy (k^2)
  if check_invariant group Poly.xx .I && check_invariant group Poly.yy .I && check_invariant group Poly.zz .I then
     IO.println "  [ALLOWED] Standard Kinetic Terms (kx^2, ky^2, kz^2 * I)"

  -- Check Linear Spin Splitting (Rashba/Dresselhaus type)
  check_linear_k_sigma group

  -- Check Altermagnetic Terms (d-wave spin splitting)
  IO.println "  Checking d-wave spin-splitting terms:"

  if check_invariant group Poly.xx .z && check_invariant group Poly.yy .z then
     -- This is loose checking
     let _ := 0

  -- We need to check (kx^2 - ky^2) * sz manually since it's a linear combination
  let check_custom (f : SPG.Core.Algebra.Vec3 → ℚ) (s : SpinComp) : Bool :=
    let test_ks : List SPG.Core.Algebra.Vec3 := [![1, 0, 0], ![0, 1, 0], ![0, 0, 1], ![1, 1, 0], ![1, 0, 1], ![0, 1, 1], ![1, 2, 3]]
    test_ks.all fun k =>
      group.all fun g =>
        let val_gk := f (act_on_k g k)
        let val_k  := f k
        if s == .I then val_gk == val_k
        else
          let s_prime := act_on_spin g s
          let coeff := project_spin s_prime s
          let is_eigen :=
              (s == .x && s_prime 1 == 0 && s_prime 2 == 0) ||
              (s == .y && s_prime 0 == 0 && s_prime 2 == 0) ||
              (s == .z && s_prime 0 == 0 && s_prime 1 == 0)
          if is_eigen then val_gk * coeff == val_k else false

  let kx2_m_ky2 (k : SPG.Core.Algebra.Vec3) : ℚ := k 0 * k 0 - k 1 * k 1
  let kx_ky     (k : SPG.Core.Algebra.Vec3) : ℚ := k 0 * k 1

  if check_custom kx2_m_ky2 .z then
    IO.println "    [ALLOWED] (kx^2 - ky^2) * σz (d-wave)"
  else
    IO.println "    [FORBIDDEN] (kx^2 - ky^2) * σz"

  if check_custom kx_ky .z then
    IO.println "    [ALLOWED] kx ky * σz (d-wave)"
  else
    IO.println "    [FORBIDDEN] kx ky * σz"

  -- Check Net Magnetization
  IO.println "  Checking Net Magnetization (M):"
  if check_custom (fun _ => 1) .z then
    IO.println "    [ALLOWED] Ferromagnetism Mz (1 * σz)"
  else
    IO.println "    [FORBIDDEN] Ferromagnetism Mz (1 * σz)"

end Demo.MagneticPhases

def main : IO Unit := do
  Demo.MagneticPhases.analyze_phase "Ferromagnet (D4h, M || z)" SPG.Material.MagneticGroup.Ferromagnet_Group_D4h_z
  Demo.MagneticPhases.analyze_phase "Antiferromagnet (D4h, PT-symmetric)" SPG.Material.MagneticGroup.Antiferromagnet_Group_PT
  Demo.MagneticPhases.analyze_phase "Altermagnet (D4h)" SPG.Material.MagneticGroup.Altermagnet_Group_D4h
