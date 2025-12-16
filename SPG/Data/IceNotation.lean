/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import Mathlib.Data.Matrix.Basic
import SPG.Core.Geometry.SpatialOps
import SPG.Core.Geometry.SpinOps

namespace SPG.Data.ICE_Notation

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps

def mk_ice_element (spatial : Matrix (Fin 3) (Fin 3) ℚ) (time_reversal : Bool) : SPGElement :=
  { spatial := spatial,
    spin := if time_reversal then spin_neg_I else spin_I }

end SPG.Data.ICE_Notation
