/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Algebra.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.Module.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Algebra.Ring.Basic
import Mathlib.Algebra.Field.Defs
import Mathlib.Algebra.Ring.Rat

namespace SPG

def magnetic_action (g : SPGElement) (v : Vec3) : Vec3 :=
  Matrix.mulVec g.spin (Matrix.mulVec g.spatial v)

def electric_action (g : SPGElement) (v : Vec3) : Vec3 :=
  Matrix.mulVec g.spatial v

end SPG
