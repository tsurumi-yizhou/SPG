/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import SPG.Physics.Hamiltonian
import SPG.Data.MagneticGroups
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace Demo.Altermagnet

open SPG
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
