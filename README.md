# FREI
**Flux Reconstruction Educational Implementation**

_by **RedHatTurtle** (Fábio Mallaco Moreira)_

_Licensed under the GNU Affero General Public License v3.0_

## Intro
This is my first attempt at writing a functional Chapel code. The code is based on a computational project from ITA's
graduate course CC-299 on high resolution methods for CFD. The code attempts to simulate a 1D shock tube experiment with
several 2nd or higher order spatial discretization methods, including the Flux Reconstruction method as described by
H.T. Huyhn from NASA Glenn.

## Compiling and running instructions

```
chpl -o frei --main-module "Control" frei.chpl
./frei --inputFile=inputExample.toml
```

# Less important shit (also work in progress)



## Coding Standards
These coding rules are meant to avoid bugs and passively improve code quality.

### Declare variables at the beginning of the block they are used

### Whenever possible use explicit typing

### Always define argument intents

### Remember to initialize variables that shouldn't default
Chapel default boolean variables to false, numbers to 0 and strings to empty. Remember to initialize variables at declaration when
these are not appropriate defaults.



## Style Standards
These style rules are meant to be only as strict as necessary to keep code easy to read, comprehend and facilitate
text editor configuration.

### Indentation
Indentation must be done with spaces, not tabs, and at an increment of 2 per level.

### Avoid unnecessary blocks {}
When issuing a single command after a loop or conditional expression use `do` instead of single line blocks.

### Use descriptive naming

Give as descriptive a name as possible, within reason. Abbreviations that would be familiar to someone outside your project with
relevant domain knowledge are OK.

### Naming conventions

Files are named after the modules they contain but start with small caps for easier autocomplete on the terminal

```
Modules    - CapitalizationAndCamelCase
Procedures - all_small_caps_and_underscores
Classes    - all_small_caps_underscores_and_c
Types      - all_small_caps_underscores_and_t
Records    - all_small_caps_underscores_and_r
Parameters - ALL_CAPS_AND_UNDERSCORE
Constants  - ALL_CAPS_AND_UNDERSCORE
Variables  - smallCapsAndCamelCase
```

```
(...) - For passing arguments and defining variable type size
[...] - For array dimensions and literals
{...} - For domain literals and code blocks
```

Default to camelCase.
People's names are **always** capitalized.
Acronyms are **always** ALL CAPS.

### Keep lines under 120 characters

### Comments
-Avoid using `\* comment *\` unless writing really long comments.
-In parameter/constant/variable declarations try to keep comments at the end of the line.
-In imperative sections add comment only lines before the code referenced.

### Use double quotes "" for strings
