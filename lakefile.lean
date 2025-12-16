import Lake
open Lake DSL

package "SPG" where
  -- add package configuration options here

lean_lib "SPG" where
  -- add library configuration options here
  -- The default root is `SPG`, which matches our file `SPG.lean`

lean_exe "demo_altermagnet" where
  root := `Demo.Altermagnet

lean_exe "demo_magnetic_phases" where
  root := `Demo.MagneticPhases

lean_exe "demo_multiferroic" where
  root := `Demo.Multiferroic

lean_exe "demo_spin_splitting" where
  root := `Demo.SpinSplitting

lean_exe "demo_coupling" where
  root := `Demo.Coupling

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"
