/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import SPG.Physics.Hamiltonian
import SPG.Physics.ResidualGroup
import SPG.Data.MagneticGroups
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

open SPG

namespace Demo.Altermagnet

open SPG
open SPG.Physics
open SPG.Physics.Hamiltonian
open SPG.Data.MagneticGroups

def ice_symbol : String := "4'22 (D4h magnetic)"

end Demo.Altermagnet

def main : IO Unit := do
  IO.println s!"Generated Group Size: {SPG.Data.MagneticGroups.Altermagnet_Group_D4h.length}"
  IO.println s!"ICE Symbol (Approx): {Demo.Altermagnet.ice_symbol}"
  IO.println "Invariant k·p Hamiltonians (basis by degree ≤ 2):"

  let group := SPG.Data.MagneticGroups.Altermagnet_Group_D4h
  let blocks := SPG.Physics.Hamiltonian.invariants_vector_by_degree_solve group 2
  for (d, hs) in blocks do
    IO.println s!"  degree {d}:"
    for h in hs do
      IO.println s!"    {SPG.Physics.Hamiltonian.ham_to_string h}"

  IO.println "Invariant complex k·p Hamiltonians (basis by degree ≤ 2):"
  let cblocks := SPG.Physics.Hamiltonian.invariants_vector_by_degree_solveC group 2
  for (d, hs) in cblocks do
    IO.println s!"  degree {d}:"
    for h in hs do
      IO.println s!"    {SPG.Physics.Hamiltonian.cham_to_string h}"

  let p_dwave : SPG.Physics.Hamiltonian.Poly :=
    (SPG.Physics.Hamiltonian.kx * SPG.Physics.Hamiltonian.kx) -
      (SPG.Physics.Hamiltonian.ky * SPG.Physics.Hamiltonian.ky)
  let p_xy : SPG.Physics.Hamiltonian.Poly :=
    SPG.Physics.Hamiltonian.kx * SPG.Physics.Hamiltonian.ky
  let p_kin : SPG.Physics.Hamiltonian.Poly :=
    (SPG.Physics.Hamiltonian.kx * SPG.Physics.Hamiltonian.kx) +
      (SPG.Physics.Hamiltonian.ky * SPG.Physics.Hamiltonian.ky)

  if SPG.Physics.Hamiltonian.isInvariantHam group (SPG.Physics.Hamiltonian.singleTerm p_dwave .z) then
    IO.println "  (kx^2 - ky^2) * σz  [d-wave altermagnetism (x^2-y^2 type)]"
  else
    IO.println "  (kx^2 - ky^2) * σz  [FORBIDDEN]"

  if SPG.Physics.Hamiltonian.isInvariantHam group (SPG.Physics.Hamiltonian.singleTerm p_xy .z) then
    IO.println "  kx ky * σz          [d-wave altermagnetism (xy type)]"
  else
    IO.println "  kx ky * σz          [FORBIDDEN]"

  if SPG.Physics.Hamiltonian.isInvariantHam group (SPG.Physics.Hamiltonian.singleTerm p_kin .I) then
    IO.println "  (kx^2 + ky^2) * I   [Standard kinetic term]"
  else
    IO.println "  (kx^2 + ky^2) * I   [FORBIDDEN]"

  -- Test residual group and Gamma point Zeeman term criterion
  IO.println "\nResidual group analysis for electric field directions:"

  let e_z : Vec3 := fun i => if i = 2 then 1 else 0  -- E || z
  let e_x : Vec3 := fun i => if i = 0 then 1 else 0  -- E || x

  IO.println "Field direction: E || z"
  let gE_z := SPG.Physics.residual_group group e_z
  IO.println s!"  Residual group order: {gE_z.length}"
  let forbids_z := SPG.Physics.forbids_gamma_zeeman_by_s group e_z
  IO.println s!"  Criterion forbids Gamma Zeeman term: {forbids_z}"

  IO.println "Field direction: E || x"
  let gE_x := SPG.Physics.residual_group group e_x
  IO.println s!"  Residual group order: {gE_x.length}"
  let forbids_x := SPG.Physics.forbids_gamma_zeeman_by_s group e_x
  IO.println s!"  Criterion forbids Gamma Zeeman term: {forbids_x}"
