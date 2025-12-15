import Lake
open Lake DSL

package "SPG" where
  -- add package configuration options here

lean_lib "SPG" where
  -- add library configuration options here
  -- The default root is `SPG`, which matches our file `SPG.lean`

lean_exe "test_basic_ops" where
  root := `Test.BasicOps

lean_exe "test_symmetry_breaking" where
  root := `Test.SymmetryBreaking

lean_exe "demo_altermagnet" where
  root := `Demo.Altermagnet

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"
