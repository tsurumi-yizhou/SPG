/-
Copyright (c) 2024 Yizhou Tong. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yizhou Tong
-/
import SPG.Core.Algebra.AlgebraBasics
import SPG.Core.Algebra.Group
import SPG.Core.Algebra.Actions
import SPG.Core.Geometry.SpatialOps
import SPG.Core.Geometry.SpinOps
import SPG.Model.SymmetryBreaking
import SPG.Model.Hamiltonian
import SPG.Model.Multiferroic
import SPG.Model.Coupling
import SPG.Model.Perturbation
import SPG.Material.MagneticGroup
import SPG.Data.IceNotation
import SPG.Data.Tetragonal
import SPG.Interface.Notation

-- Re-export common namespaces for easy usage
-- When users `import SPG`, they get access to all core modules.
