import SPG.Algebra.Basic
import Mathlib.Data.Matrix.Basic
import SPG.Geometry.SpinOps

namespace SPG.Data

open SPG.Geometry.SpinOps

def mk_ice_element (spatial : Matrix (Fin 3) (Fin 3) ℚ) (time_reversal : Bool) : SPGElement :=
  { spatial := spatial, 
    spin := if time_reversal then spin_neg_I else spin_I }

end SPG.Data
