import Mathlib.Data.Matrix.Basic

namespace SPG.Geometry.SpatialOps

def mat_4bar_z : Matrix (Fin 3) (Fin 3) ℚ :=
  ![![0, 1, 0], ![-1, 0, 0], ![0, 0, -1]]

def mat_2_x : Matrix (Fin 3) (Fin 3) ℚ :=
  ![![1, 0, 0], ![0, -1, 0], ![0, 0, -1]]

def mat_m_xy : Matrix (Fin 3) (Fin 3) ℚ :=
  ![![0, 1, 0], ![1, 0, 0], ![0, 0, 1]]

end SPG.Geometry.SpatialOps
