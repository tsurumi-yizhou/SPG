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
import SPG.Model.Coupling
import SPG.Material.MagneticGroup
import Mathlib.Data.Matrix.Basic

namespace Demo.Coupling

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps
open SPG.Interface
open SPG.Model.Coupling
open SPG.Material.MagneticGroup

-- =================================================================
-- 1. Improper Ferroelectricity (Nat. Commun. 2017)
-- =================================================================
-- "Improper electric polarization in simple perovskite oxides with two magnetic sublattices"
-- Usually involves coupling like P ~ M_A x M_B.
-- Let's check if P_z * M_Ax * M_By is allowed in a typical perovskite magnetic group.
-- Example: G-type AFM in cubic perovskite? Or lower symmetry?
-- Let's use a simple Orthorhombic group Pmmm with AFM order.
-- Or stick to D4h AFM which is common.

def run_improper_ferro_analysis : IO Unit := do
  IO.println "\n[1] Improper Ferroelectricity (Nat. Commun. 2017)"
  let group := Antiferromagnet_Group_PT -- D4h PT-symmetric
  IO.println s!"    Group: D4h (PT-symmetric AFM), Order {group.length}"

  let terms := analyze_improper_ferro group
  if terms.isEmpty then
    IO.println "    [RESULT] No trilinear coupling P M M found."
  else
    IO.println s!"    [RESULT] Found allowed couplings: {terms}"
    if terms.contains "Pz Mx My" || terms.contains "Pz My Mx" then
      IO.println "    => Pz ~ Mx * My coupling is ALLOWED."

-- =================================================================
-- 2. DM-like Interaction in Ferroelectrics (Nat. Mater. 2021)
-- =================================================================
-- "Dzyaloshinskii–Moriya-like interaction in ferroelectrics"
-- P . (L x L') or similar.
-- Symmetry-wise, this is identical to P M M coupling check.
-- We check if P (Polar) can couple to two Axial vectors antisymmetrically.
-- Actually, the DM interaction D . (S1 x S2) means Energy ~ D_k S_i S_j epsilon_ijk.
-- If D is identified with P, then Energy ~ P_k S_i S_j epsilon_ijk.
-- This requires P to transform like D (which is Axial? No, D is a vector).
-- Wait. D vector in DM interaction is usually determined by symmetry rules (Moriya rules).
-- D is an axial vector?
-- Hamiltonian H = D . (S1 x S2).
-- S1 x S2 is (Axial x Axial) -> Axial.
-- Energy must be scalar. So D . Axial = Scalar.
-- So D must be Axial? No.
-- Inversion: S->S, S->S => SxS -> S.
-- So S1xS2 is Even under I.
-- Energy is Even. So D must be Even under I.
-- D is Axial vector? Yes.
-- BUT, in the paper, they relate it to P (Polar).
-- This implies broken Inversion symmetry is key.
-- If I is broken, P exists.
-- The term is likely P . (something).
-- If "something" is Lifshitz invariant L dL, L is polar (structural)?
-- If L is polar (e.g. AFD mode), L x dL is Axial?
-- Let's assume the paper talks about P coupling to two Magnetic modes.
-- Let's just run the general coupling analysis.

def analyze_dm_like : IO Unit := do
  IO.println "\n[2] DM-like Interaction (Nat. Mater. 2021)"
  -- Let's check a polar group, e.g. C4v (Ferroelectric).
  -- Generators: C4z, Mirror_x. (No Inversion).

  let mat_4_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, -1, 0], ![1, 0, 0], ![0, 0, 1]]
  let mat_mx  : Matrix (Fin 3) (Fin 3) ℚ := ![![ -1, 0, 0], ![0, 1, 0], ![0, 0, 1]]

  let gen_C4z := Op[mat_4_z, ^1]
  let gen_Mx  := Op[mat_mx, ^1]

  let FE_Group := generate_group [gen_C4z, gen_Mx]

  IO.println s!"    Group: C4v (Ferroelectric), Order {FE_Group.length}"

  let terms := analyze_improper_ferro FE_Group -- Checking P M M
  IO.println s!"    Allowed P-M-M couplings: {terms}"

  -- The paper might refer to P coupling to Structural modes (Polar-Polar-Polar?)
  -- Let's check Polar-Polar-Polar trilinear term.
  -- e.g. P_z * u_x * u_y ?
  let check_ppp (i j k : Fin 3) := check_trilinear_invariant FE_Group .Polar .Polar .Polar i j k

  if check_ppp 2 0 1 then -- z x y
    IO.println "    [CHECK] Pz * ux * uy is ALLOWED."
  else
    IO.println "    [CHECK] Pz * ux * uy is FORBIDDEN."

end Demo.Coupling

def main : IO Unit := do
  Demo.Coupling.run_improper_ferro_analysis
  Demo.Coupling.analyze_dm_like
