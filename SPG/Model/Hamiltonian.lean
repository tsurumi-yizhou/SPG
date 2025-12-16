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
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace SPG.Model.Hamiltonian

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps
open SPG.Interface

-- General Polynomial Term: represented by a list of indices (0=x, 1=y, 2=z)
-- e.g. kx*ky^2 -> [0, 1, 1]
structure PolyTerm where
  indices : List (Fin 3)
  deriving Repr, DecidableEq, Inhabited

def mk_poly (l : List ℕ) : PolyTerm :=
  { indices := l.map (fun i => if h : i < 3 then ⟨i, h⟩ else ⟨0, by simp⟩) }

-- Common terms
def Poly.const : PolyTerm := { indices := [] }
def Poly.x : PolyTerm := mk_poly [0]
def Poly.y : PolyTerm := mk_poly [1]
def Poly.z : PolyTerm := mk_poly [2]
def Poly.xx : PolyTerm := mk_poly [0, 0]
def Poly.yy : PolyTerm := mk_poly [1, 1]
def Poly.zz : PolyTerm := mk_poly [2, 2]
def Poly.xy : PolyTerm := mk_poly [0, 1]
def Poly.yz : PolyTerm := mk_poly [1, 2]
def Poly.zx : PolyTerm := mk_poly [2, 0]
-- Cubic terms (examples)
def Poly.xxx : PolyTerm := mk_poly [0, 0, 0]
def Poly.xyz : PolyTerm := mk_poly [0, 1, 2]

def eval_poly (p : PolyTerm) (k : Vec3) : ℚ :=
  p.indices.foldl (fun acc i => acc * k i) 1

def generate_polys (degree : Nat) : List PolyTerm :=
  let rec aux (d : Nat) : List (List (Fin 3)) :=
    match d with
    | 0 => [[]]
    | d' + 1 =>
      let prev := aux d'
      prev.flatMap fun l => [l ++ [0], l ++ [1], l ++ [2]]

  -- Generate all up to degree
  let all_indices := (List.range (degree + 1)).flatMap fun d => aux d
  all_indices.map fun l => { indices := l }

-- Basis for spin matrices (sigma_x, sigma_y, sigma_z)
inductive SpinComp
| I | x | y | z
deriving Repr, DecidableEq, Inhabited

-- Action on k-vector (spatial part)
def act_on_k (g : SPGElement) (k : Vec3) : Vec3 :=
  Matrix.mulVec g.spatial k

-- Action on spin component (spin part)
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

def check_invariant (group : List SPGElement) (p : PolyTerm) (s : SpinComp) : Bool :=
  group.all fun g =>
    -- Symmetry constraint: H(k) = U H(g^-1 k) U^dagger
    -- H(k) = f(k) * sigma_s
    -- Transform: f(g^-1 k) * (g sigma_s g^-1)
    -- We check: f(g k) * (transformed sigma) == f(k) * sigma_s
    -- Note: using g k instead of g^-1 k for group average equivalence.

    let test_ks : List Vec3 := [![1, 0, 0], ![0, 1, 0], ![0, 0, 1], ![1, 1, 0], ![1, 0, 1], ![0, 1, 1], ![1, 2, 3]]

    test_ks.all fun k =>
      let val_gk := eval_poly p (act_on_k g k)
      let val_k  := eval_poly p k

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

def find_invariants (group : List SPGElement) (max_degree : Nat := 2) : List (PolyTerm × SpinComp) :=
  let terms := (generate_polys max_degree).product [.I, .x, .y, .z]
  let simple_invariants := terms.filter fun (p, s) => check_invariant group p s
  simple_invariants

-- Helper to print detailed symmetry analysis
def analyze_term_symmetry (group : List SPGElement) (f : Vec3 → ℚ) (s : SpinComp) (f_name : String) (s_name : String) : IO Unit := do
  IO.println s!"\n  Analyzing term: {f_name} * {s_name}"
  let test_k : Vec3 := ![1, 2, 3] -- Arbitrary k vector
  let val_k := f test_k

  -- Check each generator (or all elements if small)
  -- For brevity, let's just check generators + a few elements
  let elements_to_check := group.take 8 -- Check first 8 elements

  let mut i := 0
  for g in elements_to_check do
    let val_gk := f (act_on_k g test_k)

    let s_prime := act_on_spin g s
    let coeff := project_spin s_prime s

    let is_eigen :=
             (s == .x && s_prime 1 == 0 && s_prime 2 == 0) ||
             (s == .y && s_prime 0 == 0 && s_prime 2 == 0) ||
             (s == .z && s_prime 0 == 0 && s_prime 1 == 0) ||
             (s == .I)

    -- Calculate expected transformed value: val_gk * coeff
    -- It should equal val_k if invariant
    let transformed_val := val_gk * coeff
    let invariant := is_eigen && (transformed_val == val_k)

    if !invariant then
      IO.println s!"    [Broken by g{i}]"
      IO.println s!"      g{i} spatial: {repr g.spatial}"
      IO.println s!"      g{i} spin: {if g.spin == spin_neg_I then "-I (Time Reversal)" else "I"}"
      IO.println s!"      k -> g k: {repr test_k} -> {repr (act_on_k g test_k)}"
      IO.println s!"      f(k) = {val_k}, f(g k) = {val_gk}"
      IO.println s!"      σ -> g σ g⁻¹: {s_name} -> {if s == .I then "I" else s!"{repr s_prime}"}"
      IO.println s!"      Factor from spin rot: {coeff}"
      IO.println s!"      Check: {val_gk} * {coeff} ?= {val_k} => {transformed_val == val_k}"
      if !is_eigen then
        IO.println s!"      (Spin mixing occurred: {repr s_prime} is not parallel to {s_name})"

    i := i + 1

end SPG.Model.Hamiltonian
