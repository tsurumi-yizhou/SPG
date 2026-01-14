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
import SPG.Data.MagneticGroups
import SPG.Data.Tetragonal
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace Demo.MagneticPhases

open SPG
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface
open SPG.Algebra
open SPG.Physics.Hamiltonian
open SPG.Data.MagneticGroups
open SPG.Data.Tetragonal

def check_linear_k_sigma (group : List SPGElement) : IO Unit := do
  IO.println "  Checking linear spin-splitting terms (k * sigma):"
  let terms := [
    (PolyTerm.x, SpinComp.x, "kx * σx"), (PolyTerm.x, SpinComp.y, "kx * σy"), (PolyTerm.x, SpinComp.z, "kx * σz"),
    (PolyTerm.y, SpinComp.x, "ky * σx"), (PolyTerm.y, SpinComp.y, "ky * σy"), (PolyTerm.y, SpinComp.z, "ky * σz"),
    (PolyTerm.z, SpinComp.x, "kz * σx"), (PolyTerm.z, SpinComp.y, "kz * σy"), (PolyTerm.z, SpinComp.z, "kz * σz")
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

  IO.println "  Invariant k·p Hamiltonians (basis by degree ≤ 2):"
  let blocks := SPG.Physics.Hamiltonian.invariants_vector_by_degree_solve group 2
  for (d, hs) in blocks do
    IO.println s!"    degree {d}:"
    for h in hs do
      IO.println s!"      {SPG.Physics.Hamiltonian.ham_to_string h}"

  IO.println "  Invariant complex k·p Hamiltonians (basis by degree ≤ 2):"
  let cblocks := SPG.Physics.Hamiltonian.invariants_vector_by_degree_solveC group 2
  for (d, hs) in cblocks do
    IO.println s!"    degree {d}:"
    for h in hs do
      IO.println s!"      {SPG.Physics.Hamiltonian.cham_to_string h}"

  IO.println "  Invariant Hermitian complex k·p Hamiltonians (basis by degree ≤ 2):"
  let hcblocks := SPG.Physics.Hamiltonian.hermitian_invariants_vector_by_degree_solveC group 2
  for (d, hs) in hcblocks do
    IO.println s!"    degree {d}:"
    for h in hs do
      IO.println s!"      {SPG.Physics.Hamiltonian.cham_to_string h}"

  -- Check Chemical Potential
  if check_invariant group .const .I then
    IO.println "  [ALLOWED] Chemical Potential (1 * I)"
  else
    IO.println "  [FORBIDDEN] Chemical Potential (1 * I) - Strange!"

  -- Check Standard Kinetic Energy (k^2)
  if check_invariant group .xx .I && check_invariant group .yy .I && check_invariant group .zz .I then
     IO.println "  [ALLOWED] Standard Kinetic Terms (kx^2, ky^2, kz^2 * I)"

  -- Check Linear Spin Splitting (Rashba/Dresselhaus type)
  check_linear_k_sigma group

  -- Check Altermagnetic Terms (d-wave spin splitting)
  IO.println "  Checking d-wave spin-splitting terms:"
  let p_dwave : Poly := (kx * kx) - (ky * ky)
  let p_xy : Poly := kx * ky

  if isInvariantHam group (singleTerm p_dwave .z) then
    IO.println "    [ALLOWED] (kx^2 - ky^2) * σz (d-wave)"
  else
    IO.println "    [FORBIDDEN] (kx^2 - ky^2) * σz"

  if isInvariantHam group (singleTerm p_xy .z) then
    IO.println "    [ALLOWED] kx ky * σz (d-wave)"
  else
    IO.println "    [FORBIDDEN] kx ky * σz"

  -- Check Net Magnetization
  IO.println "  Checking Net Magnetization (M):"
  if isInvariantHam group (singleTerm (1 : Poly) .z) then
    IO.println "    [ALLOWED] Ferromagnetism Mz (1 * σz)"
  else
    IO.println "    [FORBIDDEN] Ferromagnetism Mz (1 * σz)"

end Demo.MagneticPhases

def main : IO Unit := do
  Demo.MagneticPhases.analyze_phase "Ferromagnet (D4h, M || z)" SPG.Data.MagneticGroups.Ferromagnet_Group_D4h_z
  Demo.MagneticPhases.analyze_phase "Antiferromagnet (D4h, PT-symmetric)" SPG.Data.MagneticGroups.Antiferromagnet_Group_PT
  Demo.MagneticPhases.analyze_phase "Altermagnet (D4h)" SPG.Data.MagneticGroups.Altermagnet_Group_D4h
  let laue := SPG.Data.Tetragonal.Laue_D2d_gens
  let mag : List SPG.SPGElement := [SPG.Data.mk_ice_element SPG.Geometry.SpatialOps.mat_4bar_z true]
  let combined := SPG.Physics.Hamiltonian.group_from_laue_magnetic laue mag
  Demo.MagneticPhases.analyze_phase "Tetragonal Laue (D2d+I) + (4bar_z*T)" combined
