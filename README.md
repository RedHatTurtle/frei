# FREI
**Flux Reconstruction Educational Implementation**

_by **Fábio Malacco Moreira** (RedHatTurtle)_

_Licensed under the GNU Affero General Public License v3.0_

## Intro
This code is a prototype for a high-order Flux Reconstruction method (developed by H.T. Huyhn from NASA Glenn) LES
solver. Secondarily, it should also be able to incorporate RANS modeling and a Finite Volume solver.

The initial objective is to build an accurate, stable and efficient 3D Euler solver while also providing the basic
building blocks and structure for future development of the solver.

Currently the code has basic support for 1D and 2D inviscid, smooth flows with an experimental 3D implementation under
development. An experimental limiter for transonic flows is implemented but not maintained.

The next major objectives are implementing multi locale parallelism and benchmarking/tuning of the code up to acceptable
performance standards compared to similar solver in Fortran/C/C++.

## Compiling and running instructions

```
chpl -o frei --main-module "Control" frei.chpl
./frei --inputFile=inputExample.toml
```

# Less important shit (also work in progress)

## Git Standards

### Branch Naming

Main branches:
- `main`: Main trunk, most ready for users versions of the code, Stable.
- `dev` : Current development version, Unstable.

Default branch naming: `<author>/<branch-type>/<issueID>-<branch-name>`

- `Author`: Short author/owner identifier. Ex: `rht`, `fabio`, `fmoreira`.
- `Issue ID`: ID# of bug fix, feature request, etc, if existing.
- `Branch Name`: Short descriptive name of work done, ~3 words in camelCase.
- `Branch Type`: One of
  + `fix`: Bug fixes on main trunk.
  + `wip`: Work in Progress, for unfinished development branches.
  + `feat`: Finished new feature development branch, to be or already merged.
  + `dbg`: Debug versions of the code, never expected to be merged to main.
  + `test`: Test case configurations update only.
  + `patch`: Updates for compatibility reasons with no intended behavior change.
  + `docs` : Update to documentation only

style: (formatting, missing semi colons, etc; no production code change)
refactor: (refactoring production code, eg. renaming a variable)
test: (adding missing tests, refactoring tests; no production code change)

Examples:
- `rht/fix/47-isothermalWallEnergy` Fixing the energy gradient on an isothermal wall boundary condition
- `fabio/wip/compressionLimiter` Developing a new limiter
- `fmoreira/patch/chapel-1.30.0` Updating code for new compiler version

### Merge Strategy

#### No Fast-Forward Merges
Disabling fast-forward merges to the `main` and `dev` allows a cleaner history for these branches while keeping the development history of the merged code in it's original branch.

## Chapel Standards
These coding rules are meant to avoid bugs and passively improve code quality and readability.

### Declare variables at the beginning of the code block they are used

### Whenever possible use explicit typing

### Always define argument intents

### Remember to initialize variables that shouldn't start with the type's default value
Chapel default boolean variables to false, numbers to 0 and strings to empty. Remember to initialize variables at
declaration when these are not appropriate defaults.

## Style Standards
These style rules are meant to be only as strict as necessary to keep code easy to read, comprehend and facilitate text
editor configuration.

### Indentation
Indentation must be done with spaces, not tabs, and at an increment of 2 per level.

### Avoid unnecessary blocks {}
When issuing a single command after a loop or conditional expression favor the single statement form of these
expressions (ex: `for ... do`, `if ... then`) instead of blocks with only one statement.

### Use descriptive naming
Give as descriptive a name as possible, within reason. Abbreviations that would be familiar to someone outside your
project with relevant domain knowledge are OK.

### Naming conventions
Files are named after the modules they contain but start with small caps for easier auto-complete on the terminal

```
Modules    - PascalCase
Parameters - SCREAMING_SNAKE_CASE
Constants  - SCREAMING_SNAKE_CASE
Variables  - camelCase
Domains    - sameNameAsArrayAnd_d
Procedures - snake_case
Classes    - snake_case_and_c
Types      - snake_case_and_t
Records    - snake_case_and_r
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
- Avoid using `\* comment *\` unless writing really long comments or required by **ChplDoc**.
- In parameter/constant/variable declarations try to keep comments at the end of the line.
- In imperative sections add comment only lines before the code referenced.

### Use double quotes "" for strings

## Order vs Degree

When referring to the solution interpolation polynomial it is usual in the High-Order CFD community to use the term
"order" to refer to the degree of the polynomial. This can lead to miscommunication and misunderstanding in discussions,
as many of my colleagues can attest. This is in no way a topic of great concern currently but it is a definition yet to
be made. Currently my predisposition is to use the term "degree", leaving the term "order" only for the few places where
the "order of approximation", which is related to the approximation error, is being referred to.

Currently the terms "order" and "degree" may be used interchangeably in the code and are always (mostly?) referring to
the degree of the polynomial. Ex: `iOrder` (interpolation order), `solOrder` (solution order).

References:
- [Wikipedia entry for "degree of a polynomial"](https://en.wikipedia.org/wiki/Degree_of_a_polynomial)
- [Wikipedia entry for "order of approximation"](https://en.wikipedia.org/wiki/Order_of_approximation)

## Common Abbreviations
It seems like a good rule of thumb to try to use 4 letter abbreviations.

### Mesh related terms
- `node`, `line`, `tria`, `quad`, `tetr`, `pyra`, `pris`, `hexa` : Mesh element topologies/geometries
- `node`, `edge`, `face`, `cell` : Mesh Elements Hierarchy
- `vert` : Vertices, nodes that are in the corner of an edge/face/cell
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

- `xStat` : Static X
- `xDyna` : Dynamic X
- `xStag` : Stagnation X

- `Dens`   : Density
- `Ener`   : Energy
- `Pres`   : Static Pressure
- `Temp`   : Temperature
- `VelV`, `VelX`, `VelY`, `VelZ` : Velocity vector and it's components
- `MomV`, `MomX`, `MomY`, `MomZ` : Momentum vector and it's components
- `VelM, VelMag        `   : Velocity magnitude
- `Mach`   : Mach Number
- `VelA`, `A` : Speed of Sound
- `EneK`, `EnerKin`, `ek` : Kinetic Energy / Dynamic Pressure
- `EneI`, `EnerInt`, `ei` : Internal Energy
- `Enth`, `H` : Enthalpy
- `Entr`, `S` : Entropy

- `varFree` : Free stream property
- `varRef`  : Reference value for a property, usually the free stream conditions
- `varScal` : Nondimensionalization scale

### FR Specific terms
- `Sol`, `Con`, `Cons`, `cv` : Conserved variables
- `Flx`,        `Flux`, `fv` : Flux vector
- `Pri`,        `Prim`, `pv` : Primitive variables

- `SolPnt`, `SP` : Solution Point
- `FlxPnt`, `FP` : Flux Point
- `Point` , `PT` : Generic, possibly auxiliary, interpolation point

- `meshSP`, `meshFP` : Global mesh index of a mesh element (SP, FP, etc)
- `faceFP`, `faceNode` : Local cell index/position of a mesh element (SP, FP, Node, etc)
- `cellSP`, `cellFP` : Local cell index/position of a mesh element (SP, FP, etc)

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
- `alpha` : Rotation on the Z axis (wingspan axis, oriented left to right), aka: Angle of Attack ($\alpha$), Pitch ($\theta$), Flight path ($\gamma$).
- `beta`  : Rotation on the Y axis (vertical axis, oriented down to up), aka: Side-slip ($\beta$), Yaw ($\psi$), Heading ($\sigma$).
- `phi`   : Rotation on the X axis (fuselage axis, oriented front to back), aka: Roll ($\phi$), Bank ($\mu$).

- `lift` : Lift ($L$)
- `drag` : Drag ($D$)
- `latf` : Lateral force ($Q$)

## References

Input files are written in [TOML](https://toml.io/en/)
