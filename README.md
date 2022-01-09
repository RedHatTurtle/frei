# FREI
**Flux Reconstruction Educational Implementation**

_by **Fábio Malacco Moreira** (RedHatTurtle)_

_Licensed under the GNU Affero General Public License v3.0_

## Intro
This code is a prototype version of a full high-order Flux Reconstruction (developed by H.T. Huyhn from NASA Glenn) LES
solver.

The initial objective is to build an accurate, stable and efficient 1D Euler solver while also providing the basic
building blocks and structure for future development of the solver. Most of the work done currently is still targeted
at this goal. The next major target will be to expand the solver to two dimensional flows and add support for
simulating laminar viscous flows described by the Navier-Stokes equations.

## Compiling and running instructions

```
chpl -o frei --main-module "Control" frei.chpl
./frei --inputFile=inputExample.toml
```

# Less important shit (also work in progress)

## Coding Standards
These coding rules are meant to avoid bugs and passively improve code quality and readability.

### Declare variables at the beginning of the code block they are used

### Whenever possible use explicit typing

### Always define argument intents

### Remember to initialize variables that shouldn't start with the type's default value
Chapel default boolean variables to false, numbers to 0 and strings to empty. Remember to initialize variables at
declaration when these are not appropriate defaults.

## Style Standards
These style rules are meant to be only as strict as necessary to keep code easy to read, comprehend and facilitate
text editor configuration.

### Indentation
Indentation must be done with spaces, not tabs, and at an increment of 2 per level.

### Avoid unnecessary blocks {}
When issuing a single command after a loop or conditional expression use `do` instead of single line blocks.

### Use descriptive naming
Give as descriptive a name as possible, within reason. Abbreviations that would be familiar to someone outside your
project with relevant domain knowledge are OK.

### Naming conventions
Files are named after the modules they contain but start with small caps for easier auto-complete on the terminal

```
Modules    - CapitalizationAndCamelCase
Parameters - ALL_CAPS_AND_UNDERSCORE
Constants  - ALL_CAPS_AND_UNDERSCORE
Variables  - smallCapsAndCamelCase
Domains    - sameNameAsArrayAnd_d
Procedures - all_small_caps_and_underscores
Classes    - all_small_caps_underscores_and_c
Types      - all_small_caps_underscores_and_t
Records    - all_small_caps_underscores_and_r
```

```
(...) - For passing arguments and defining variable type size
[...] - For array dimensions and literals
{...} - For domain literals and code blocks
```

Default to camelCase.
People's names are **always** Capitalized.
Acronyms are **always** ALL CAPS.

### Limit code lines to 120 columns

### Limit comment/text lines to 120 characters

### Comments
-Avoid using `\* comment *\` unless writing really long comments.
-In parameter/constant/variable declarations try to keep comments at the end of the line.
-In imperative sections add comment only lines before the code referenced.

### Use double quotes "" for strings

## Order vs Degree

When referring to the solution interpolation polynomial it is usual in the High-Order CFD community to use the term
"order" to refer to the degree of the polynomial. This can lead to miscommunication and misunderstanding in
discussions, as many of my colleagues can attest. This is in no way a topic of great concern currently but it is a
definition yet to be made. Currently my predisposition is to use the term "degree", leaving the term "order" only for
the few places where the "order of approximation", which is related to the approximation error, is being referred to.

Currently the terms "order" and "degree" may be used interchangeably in the code and are always (mostly?) referring to
the degree of the polynomial. Ex: `iOrder` (interpolation order), `solOrder` (solution order).

References:
- [Wikipedia entry for "degree of a polynomial"](https://en.wikipedia.org/wiki/Degree_of_a_polynomial)
- [Wikipedia entry for "order of approximation"](https://en.wikipedia.org/wiki/Order_of_approximation)

## Common Abbreviations
It seems like a good rule of thumb to try to use 4 letter abbreviations.

### Mesh related terms
- `node`, `line`, `tria`, `quad`, `tetr`, `pyra`, `pris`, `hexa` : Mesh element topologies/geometries
- `node`, `edge`, `face`, `cell` : Mesh Elements
- `boco`, `bc` : Boundary Condition
- `faml` : Family
- `ghst` : Ghost cell/FP
- `host` : The cell neighboring a ghost

- `uni` : Unit vector
- `nrm` : Face normal vector

### Physical model / equation set related terms
- `Invs` : Inviscid
- `Visc` : Viscous

- `Sub` : Subsonic
- `Sup` : Supersonic

- `Dens`   : Density
- `Ener`   : Energy
- `Pres`   : Static Pressure
- `Temp`   : Temperature
- `Vel`    : Velocity
- `Mom`    : Momentum
- `Mach`   : Mach Number
- `A`      : Speed of Sound
- `EneK`, `EneKin`, `ek` : Kinetic Energy / Dynamic Pressure
- `EneI`, `EneInt`, `ei` : Internal Energy
- `Enth`, `H` : Enthalpy

### FR Specific terms
- `Sol`, `Con`, `Cons`, `cv` : Conserved variables
- `Flx`,        `Flux`, `fv` : Flux vector
- `Pri`,        `Prim`, `pv` : Primitive variables

- `SolPnt`, `SP` : Solution Point
- `FlxPnt`, `FP` : Flux Point
- `Point` , `PT` : Generic interpolation point

- `Phys` : Physical domain
- `Comp` : Computational domain

- `xyz[]`= $[ x ,  y  ,  z   ]$ : Vector of physical coordinates
- `rst[]`= $[\xi, \eta, \zeta]$ : Vector of computational coordinates

- `Eq`  : Equation
- `Idx` : Index
- `Cnt` : Count

### Mathematical operators

- `Interp` : Interpolation
- `Deriv`  : Derivative
- `Grad`   : Gradient
- `Div`    : Divergence
- `Rot`    : Rotational

### GFR Boundary Condition functions abbreviations
- `gufp` : Ghost solution (U) at Flux Point
- `hufp` : Host solution (U) at Flux Point
- `fxyz` : Coordinates of the Flux Points
- `fnrm` : Face Normal
- `cv_in` : Conserved Variables
- `pv_in` : Primitive Variables
- `wall_temperature` : Wall Temperature
- `var_int` : Variable in the **interior FP** (host cell)
- `var_ext` : Variable in the **exterior FP** (ghost cell)
- `var_ref` : Variable prescribed for the boundary condition
- `var0`    : Total var. Ex: Total Pressure `p0`, Total Temperature `t0`

### Flight Mechanics
- `alpha` : Rotation on the Y axis (wingspan axis), aka: Angle of Attack ($\alpha$), Pitch ($\theta$), Flight path ($\gamma$).
- `beta`  : Rotation on the Z axis (vertical axis), aka: Sideslip ($\beta$), Yaw ($\psi$), Heading ($\sigma$).
- `phi`   : Rotation on the X axis (fuselage axis), aka: Roll ($\phi$), Bank ($\mu$).

- `lift` : Lift ($L$)
- `drag` : Drag ($D$)
- `latf` : Lateral force ($Q$)

## References

Input files are written in [TOML](https://toml.io/en/)
