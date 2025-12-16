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

namespace SPG.Material.MagneticGroup

open SPG.Core.Algebra
open SPG.Core.Geometry.SpatialOps
open SPG.Core.Geometry.SpinOps
open SPG.Interface

-- Define common matrices
def mat_id : Matrix (Fin 3) (Fin 3) ℚ := 1
def mat_inv : Matrix (Fin 3) (Fin 3) ℚ := -1
def mat_4_z : Matrix (Fin 3) (Fin 3) ℚ := ![![0, -1, 0], ![1, 0, 0], ![0, 0, 1]]
def mat_2_x : Matrix (Fin 3) (Fin 3) ℚ := ![![1, 0, 0], ![0, -1, 0], ![0, 0, -1]]
def mat_2_xy : Matrix (Fin 3) (Fin 3) ℚ := ![![0, 1, 0], ![1, 0, 0], ![0, 0, -1]]

-- 1. Ferromagnet (FM)
def gen_FM_C4z : SPGElement := Op[mat_4_z, ^1]
def gen_FM_Inv : SPGElement := Op[mat_inv, ^1]

def Ferromagnet_Group_D4h_z : List SPGElement := generate_group [gen_FM_C4z, gen_FM_Inv]


-- 2. Antiferromagnet (AFM)
def gen_AFM_C4z : SPGElement := Op[mat_4_z, ^1]
def gen_AFM_PT  : SPGElement := Op[mat_inv, ^-1] -- P * T
def gen_AFM_C2x : SPGElement := Op[mat_2_x, ^1] -- Rotation around x

def Antiferromagnet_Group_PT : List SPGElement := generate_group [gen_AFM_C4z, gen_AFM_PT, gen_AFM_C2x]


-- 3. Altermagnet (AM) - D4h (from Demo)
def gen_AM_C4z_TR : SPGElement := Op[mat_4_z, ^-1]
def gen_AM_C2xy   : SPGElement := Op[mat_2_xy, ^1]
def gen_AM_Inv    : SPGElement := Op[mat_inv, ^1]

def Altermagnet_Group_D4h : List SPGElement := generate_group [gen_AM_C4z_TR, gen_AM_C2xy, gen_AM_Inv]

end SPG.Material.MagneticGroup
