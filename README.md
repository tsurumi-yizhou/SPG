# SPG (Spin Point Groups) Library

SPG is a Lean 4 library for analyzing Spin Point Groups, Magnetic Point Groups (MPGs), and symmetry breaking phenomena in condensed matter physics. It provides tools to define magnetic symmetries, generate group closures, and verify physical properties such as allowed electric polarization.

## Installation

Add the following dependency to your `lakefile.lean`:

```lean
require SPG from git
  "https://github.com/tsurumi-yizhou/SPG.git"
```

Then update dependencies:

```bash
lake update
```

## Usage Guide

This section demonstrates how to use the library to analyze magnetic symmetries.

### 1. Import Modules

Import the main library to access all core functionalities:

```lean
import SPG

open SPG
open SPG.Interface
open SPG.Algebra
open SPG.Physics
open SPG.Geometry.SpatialOps
```

### 2. Define Group Generators

Use the `Op[spatial, spin]` syntax to define group elements.
* `^1`: Identity spin operation (non-magnetic).
* `^-1`: Time-reversal operation (spin flip).

```lean
-- Example: 4-fold rotoinversion along z (non-magnetic)
def g1 : SPGElement := Op[mat_4bar_z, ^1]

-- Example: 2-fold rotation along x combined with time-reversal
def g2 : SPGElement := Op[mat_2_x, ^-1]
```

### 3. Generate Group Closure

Use `generate_group` to compute the full group from a set of generators.

```lean
def my_group : List SPGElement := generate_group [g1, g2]

-- Check group order
#eval my_group.length
```

### 4. Symmetry Breaking Analysis

Given a magnetic ordering vector (Neel vector), compute the Magnetic Point Group (MPG). The MPG consists of all elements in the original group that leave the magnetic vector invariant.

```lean
def neel_vector : Vec3 := ![1, 1, 0]

-- Compute the Magnetic Point Group
def my_mpg : List SPGElement := get_mpg my_group neel_vector

-- Check the size of the MPG
#eval my_mpg.length
```

### 5. Physical Property Verification

Verify if the resulting symmetry group allows specific physical properties, such as spontaneous polarization along the z-axis.

```lean
-- Returns true if z-polarization is allowed, false otherwise
#eval allows_z_polarization my_mpg
```

### 6. Formal Proof

Physical conclusions can be formally verified using `native_decide`.

```lean
theorem polarization_allowed : 
  (allows_z_polarization my_mpg) = true := by
  native_decide
```

## Core Modules

* **SPG.Interface.Notation**: Syntax definitions for `Op[...]`, `^1`, etc.
* **SPG.Algebra.Group**: Algorithms for finite group generation.
* **SPG.Physics.SymmetryBreaking**: Functions for MPG extraction and property verification.
* **SPG.Geometry.SpatialOps**: Predefined spatial operation matrices (e.g., `mat_4bar_z`).
* **SPG.Demo.Altermagnet**: Example implementation for d-wave altermagnet.

## Notes

* **Precision**: All matrix and vector operations use `Rational` numbers (`ℚ`) to ensure exact results and decidability.
* **Scope**: The `generate_group` function is optimized for finite Point Groups.

## Acknowledgments

Special thanks to **Gemini 3** for the contributions to the code structure and implementation details.
