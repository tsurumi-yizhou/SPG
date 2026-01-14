# SPG (Spin Point Groups) Library

> [中文文档 (Chinese Documentation)](README_zh.md)

SPG is a Lean 4 library designed for exact algebraic analysis of magnetic point groups and their representations in condensed matter physics.

Unlike numerical tools (Python/MATLAB), SPG operates entirely over the field of rational numbers (`ℚ`), ensuring that all symmetry checks, invariant generations, and "forbidden term" classifications are mathematically exact and decidable.

## Key Features

*   **Exact Algebra**: No floating-point errors. Matrices and coefficients are `ℚ`.
*   **Invariant Analysis**: Automated solving of linear equations to find allowed k·p Hamiltonian terms (e.g., checking for altermagnetism).
*   **Symmetry Breaking**: Calculate Magnetic Point Groups (MPG) from a given magnetic order parameter.
*   **Physical Properties**: Verify macroscopic properties like allowed electric polarization or net magnetization.

## Quick Start

### Installation

Add to your `lakefile.lean`:

```lean
require SPG from git
  "https://github.com/tsurumi-yizhou/SPG.git"
```

### Example: Checking Altermagnetism

```lean
import SPG

open SPG
open SPG.Interface
open SPG.Algebra
open SPG.Physics.Hamiltonian
open SPG.Geometry.SpatialOps

-- Define a group generator: C4z * TimeReversal
def g_C4z_T : SPGElement := Op[mat_4_z, ^-1]
def g_C2x   : SPGElement := Op[mat_2_x, ^1]
def g_Inv   : SPGElement := Op[mat_inv, ^1]

-- Generate the full group
def my_group : List SPGElement := generate_group [g_C4z_T, g_C2x, g_Inv]

-- Check if d-wave spin splitting (kx^2 - ky^2)σz is allowed
-- This is the hallmark of d-wave altermagnetism
def p_dwave : Poly := (kx * kx) - (ky * ky)
#eval isInvariantHam my_group (singleTerm p_dwave .z)
```

## Documentation

For a detailed guide on the physical model, implementation details, and advanced usage, please refer to the **[Chinese Documentation](README_zh.md)** (currently the most comprehensive guide).

## Source Code Map

*   **[SPG/Algebra/Basic.lean](SPG/Algebra/Basic.lean)**: `SPGElement` definition.
*   **[SPG/Algebra/Actions.lean](SPG/Algebra/Actions.lean)**: Physical actions on polar/axial vectors.
*   **[SPG/Physics/Hamiltonian/Invariants.lean](SPG/Physics/Hamiltonian/Invariants.lean)**: Algorithms for finding invariant Hamiltonians.
