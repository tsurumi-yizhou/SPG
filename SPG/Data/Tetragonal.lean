/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Geometry.SpatialOps
import SPG.Interface.Notation
import SPG.Data.IceNotation

namespace SPG.Data.Tetragonal

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Data.ICE_Notation

-- Generators for CuFeS2 (Chalcopyrite structure, I-42d, #122)
-- Using 4bar_z, 2_x, and m_xy (standard generators for D2d point group)
-- Assuming non-magnetic for now (time_reversal = false)

def gen1 : SPGElement := mk_ice_element mat_4bar_z false
def gen2 : SPGElement := mk_ice_element mat_2_x false
def gen3 : SPGElement := mk_ice_element mat_m_xy false

end SPG.Data.Tetragonal
