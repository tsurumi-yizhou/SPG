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
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace Demo.Altermagnet

open SPG
open SPG.Geometry.SpatialOps
open SPG.Geometry.SpinOps
open SPG.Interface
open SPG.Algebra

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

-- 2. Hamiltonian Analysis
-- We want to find H(k) ~ c_ij k_i k_j terms allowed by symmetry.
-- General form: H(k) = sum_{i,j} A_{ij} k_i k_j
-- Symmetry constraint: g H(k) g^{-1} = H(g k)
-- For spin-independent hopping (kinetic energy), H is a scalar in spin space (Identity).
-- But altermagnetism involves spin-dependent terms.
-- Altermagnet Hamiltonian usually looks like: H = (k^2 terms) * I_spin + (k-dependent field) * sigma
-- Let's analyze terms of form: k_a k_b * sigma_c

-- Basis for quadratic k terms (k_x^2, k_y^2, k_z^2, k_x k_y, k_y k_z, k_z k_x)
inductive QuadTerm
| xx | yy | zz | xy | yz | zx
deriving Repr, DecidableEq, Inhabited

def eval_quad (q : QuadTerm) (k : Vec3) : ℚ :=
  match q with
  | .xx => k 0 * k 0
  | .yy => k 1 * k 1
  | .zz => k 2 * k 2
  | .xy => k 0 * k 1
  | .yz => k 1 * k 2
  | .zx => k 2 * k 0

def all_quads : List QuadTerm := [.xx, .yy, .zz, .xy, .yz, .zx]

-- Basis for spin matrices (sigma_x, sigma_y, sigma_z)
inductive SpinComp
| I | x | y | z
deriving Repr, DecidableEq, Inhabited

-- Action on k-vector (spatial part)
def act_on_k (g : SPGElement) (k : Vec3) : Vec3 :=
  Matrix.mulVec g.spatial k

-- Action on spin component (spin part)
-- sigma_prime = U sigma U^dagger
-- Since our spin matrices are +/- I, the conjugation is trivial?
-- WAIT. The user wants "d-wave altermagnet".
-- In altermagnets, the spin symmetry is usually non-trivial (e.g. spin flip).
-- But our SPGElement definition only supports spin matrices as "numbers" (Matrix 3x3).
-- Currently `spin` is just +/- I.
-- If spin part is just +/- I (Time Reversal), then:
-- T sigma T^{-1} = - sigma (Time reversal flips spin)
-- So if g.spin = -I (Time Reversal), the spin vector flips.
-- If g.spin = I, spin vector stays.

def act_on_spin (g : SPGElement) (s : SpinComp) : Vec3 :=
  let s_vec : Vec3 := match s with
    | .x => ![1, 0, 0]
    | .y => ![0, 1, 0]
    | .z => ![0, 0, 1]
    | .I => ![0, 0, 0] -- Handle I separately

  if s == .I then ![0, 0, 0] -- I transforms to I (scalar)
  else
    -- Apply spatial rotation R to the axial vector sigma
    -- sigma' = (det R) * R * sigma
    -- If spin part has T (-I), then sigma' = - sigma'
    let detR := Matrix.det g.spatial
    let rotated := Matrix.mulVec g.spatial s_vec
    let axial_rotated := detR • rotated -- Scalar multiplication

    if g.spin == spin_neg_I then
      -axial_rotated
    else
      axial_rotated

-- Helper to check if a transformed spin vector matches a target basis component (with sign)
-- Returns 0 if orthogonal, 1 or -1 if parallel/antiparallel
def project_spin (v : Vec3) (target : SpinComp) : ℚ :=
  match target with
  | .x => v 0
  | .y => v 1
  | .z => v 2
  | .I => 0

def check_invariant (q : QuadTerm) (s : SpinComp) : Bool :=
  Altermagnet_Group.all fun g =>
    -- Symmetry constraint: H(k) = U H(g^-1 k) U^dagger
    -- H(k) = f(k) * sigma_s
    -- Transform: f(g^-1 k) * (g sigma_s g^-1)
    -- We check: f(g k) * (transformed sigma) == f(k) * sigma_s
    -- Note: using g k instead of g^-1 k for group average equivalence.

    let test_ks : List Vec3 := [![1, 0, 0], ![0, 1, 0], ![0, 0, 1], ![1, 1, 0], ![1, 0, 1], ![0, 1, 1], ![1, 2, 3]]

    test_ks.all fun k =>
      let val_gk := eval_quad q (act_on_k g k)
      let val_k  := eval_quad q k

      if s == .I then
         -- Scalar term: val_gk == val_k
         val_gk == val_k
      else
        -- Vector term: val_gk * (transformed sigma) == val_k * sigma_s
        -- We only support cases where transformed sigma is parallel to sigma_s
        -- (i.e. no mixing like x -> y).
        -- Let's project transformed sigma onto s.
        let s_prime := act_on_spin g s
        let coeff := project_spin s_prime s

        -- Check if it stays in the same direction (possibly with sign change)
        -- And check if orthogonal components are zero (no mixing)
        let is_eigen :=
             (s == .x && s_prime 1 == 0 && s_prime 2 == 0) ||
             (s == .y && s_prime 0 == 0 && s_prime 2 == 0) ||
             (s == .z && s_prime 0 == 0 && s_prime 1 == 0)

        if is_eigen then
          val_gk * coeff == val_k
        else
          false -- If mixing occurs, this single term is not an invariant by itself.

def find_invariants : List (QuadTerm × SpinComp) :=
  let terms := (all_quads.product [.I, .x, .y, .z])
  -- Also consider symmetric combinations like (kx^2 + ky^2) and antisymmetric (kx^2 - ky^2)
  -- Currently we only check pure basis terms.
  -- But (kx^2 + ky^2) * I should be invariant.
  -- Let's manually check (kx^2 + ky^2) * I and (kx^2 - ky^2) * sigma_z?
  -- Or better, construct a linear combination solver?
  -- For now, let's just stick to pure basis.
  -- If kx^2 and ky^2 transform into each other, neither is invariant alone.
  -- But their sum might be.
  -- To properly find Hamiltonian, we need representation theory (projectors).
  -- Simplification: Check "Is kx^2 + ky^2 invariant?" explicitly.

  let simple_invariants := terms.filter fun (q, s) => check_invariant q s
  simple_invariants

-- 3. ICE Symbol Output (Simplistic)
-- We just print generators for now, as full ICE symbol logic is complex.
def ice_symbol : String :=
  "4'22 (D4h magnetic)" -- Placeholder, deriving from generators

end Demo.Altermagnet

def main : IO Unit := do
  IO.println s!"Generated Group Size: {Demo.Altermagnet.Altermagnet_Group.length}"
  IO.println s!"ICE Symbol (Approx): {Demo.Altermagnet.ice_symbol}"
  IO.println "Allowed Quadratic Hamiltonian Terms H(k):"

  let invariants := Demo.Altermagnet.find_invariants
  for (q, s) in invariants do
    let q_str := match q with
      | .xx => "kx^2" | .yy => "ky^2" | .zz => "kz^2"
      | .xy => "kx ky" | .yz => "ky kz" | .zx => "kz kx"
    let s_str := match s with
      | .I => "I" | .x => "σx" | .y => "σy" | .z => "σz"
    IO.println s!"  {q_str} * {s_str}"

  -- Manually check d-wave altermagnet term: (kx^2 - ky^2) * sigma_z?
  -- Or kx ky * sigma_z ?
  -- Let's define a custom checker for kx^2 - ky^2
  let check_custom (f : SPG.Vec3 → ℚ) (s : Demo.Altermagnet.SpinComp) : Bool :=
    let test_ks : List SPG.Vec3 := [![1, 0, 0], ![0, 1, 0], ![0, 0, 1], ![1, 1, 0], ![1, 0, 1], ![0, 1, 1], ![1, 2, 3]]
    test_ks.all fun k =>
      -- Re-implement loop
      Demo.Altermagnet.Altermagnet_Group.all fun g =>
        let val_gk := f (Demo.Altermagnet.act_on_k g k)
        let val_k  := f k
        -- Re-use projection logic from check_invariant
        if s == .I then
          val_gk == val_k
        else
          let s_prime := Demo.Altermagnet.act_on_spin g s
          let coeff := Demo.Altermagnet.project_spin s_prime s

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

  if check_custom kx_ky .z then
    IO.println "  kx ky * σz          [d-wave altermagnetism (xy type)]"

  if check_custom kx2_plus_ky2 .I then
    IO.println "  (kx^2 + ky^2) * I   [Standard kinetic term]"
