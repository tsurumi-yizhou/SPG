import SPG.Data.ICE_Notation
import SPG.Geometry.SpatialOps

namespace SPG.Data.Tetragonal

open SPG.Geometry.SpatialOps

-- Generators for CuFeS2 (Chalcopyrite structure, I-42d, #122)
-- Using 4bar_z, 2_x, and m_xy (standard generators for D2d point group)
-- Assuming non-magnetic for now (time_reversal = false)

def gen1 : SPGElement := mk_ice_element mat_4bar_z false
def gen2 : SPGElement := mk_ice_element mat_2_x false
def gen3 : SPGElement := mk_ice_element mat_m_xy false

end SPG.Data.Tetragonal
