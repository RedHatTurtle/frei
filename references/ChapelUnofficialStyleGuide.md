# The Unofficial Chapel Style Guide

Chapel styles still have a lot of grey areas. This guide tries to capture some of the more black and white stylistic choices Chapel users typically choose.

# Style Guide

## Casing

Many symbols are named with *camelCase* and *PascalCase*. The *under_score* and *lowercase* patterns are not preferred.

## Variables

Variable names generally follow *camelCase* naming, with the exception of single-character identifiers.

```chpl
  var foo = 1;
  var fooBar = 2;
  var A = [1,2,3];
```

Sync variable names end with a '$' to remind users of the communication costs
    associated with reading/writing the variable.
    
 ```chpl
  var sy$: sync int = 1
```

## Object-oriented programming

Classes and Records generally followed *PascalCase* naming:

```chpl
record DataPoint {
  var x, y, z: real;
}

class DataPoint {
  var x, y, z: real;
}
```

## Modules

* Modules and module file names follow *PascalCase*:

```chpl
// AwesomeModule.chpl

module AwesomeModule { }
```

## Whitespace

* Use 2 space indentation.
  * Use soft tabs set to 2 spaces. 
  * Hard tabs are frowned upon, due to the inconsistency of how they are displayed in editors.

* Avoid having lines of code longer than 80 characters (including whitespace)
    * This ensures readability and maintainability

* Leave a space after opening a comment:

    ```chapel
    /* The best variable */
    var foo = 1;
    ```

* Trailing white space is a sin. It's a good idea to make it visible in your editor, so that it's easy to catch/remove.

## Comments

Use `chpldoc`-style comments whenever documenting a declaration:

```chpl
  /* Module documentation */
  module FooBar {

  /* Variable documentation */
  var foo = 1;

  /* Procedure documentation */
  proc getFoo() { return foo; }

  /* Record documentation */
  record Bar { }

  /* Even use this style for non-publicly documented symbols */
  private proc _incrementFoo() { foo += 1; }
}
```


## Variables

Variables use *camelCase* casing, except for single character identifiers, which can be lower or upper case, depending on the context.

```chpl
var dom = {1..10, 1..10};
var distDom = dom dmapped Block(dom); 

// Capitalized, because it's a single character.
var A: [dom] int;

// Sometimes lowercase is used for single character identifier to distinguish 1D arrays (vectors) from 2D arrays (matrices)
var v: [1..10] int;
```

### Variable declarations

#### Repetitive declarations

Multiple declarations of the same type can be declared on the same line or with a multi-line declaration if desired:

Same line:
```chpl
var foo, bar, baz: int;
```

Multi-line declaration:
```chpl
var foo: int,
    bar: int,
    baz: int;
```

This style is also OK:
```chpl
var foo: int;
var bar: int;
var baz: int;
```

#### Types vs. Initial values


There is no strong preference between declaring variables with explicit values or explicit types:


Explicit values are OK:
```chpl
var x = 0;
var y = 0.0;
```

Explicit types are also OK:
```chpl
var x: int;
var y: real;
```


## Blocks

There are several patterns adopted in writing block statements.

### Braces

Brace style can vary. Try to be consistent with the file(s)/project you are working in:

* Option 1: Compact braces trades readability for saved vertical space:

    ```chapel

    if foo > bar {
      doThis();
    } else {
      doThat();
    }
    ```

* Option 2: Loose braces (with newlines) trades vertical space for improved clarity and readability:

    ```chapel

    if foo > bar
    {
      doThis();
    }
    else
    {
      doThat();
    }
    ```
 
 ### Braceless
 
 Dropping braces for a `then` or `do` keyword reduces syntax pollution and improves readability, but comes at a cost of flexibility, e.g. if we want to add a `writeln()` inside the block, we'll need to switch back to braces. 
 
 Braceless style can be used at the author's discretion.
    
* Option 1: `then\n` style:

    ```chapel
    if foo > bar then
      doThis();
    else
      doThat();
    ```
    
* Option 2: `\nthen` style:

    ```chapel
    if foo > bar
      the doThis();
    else
      doThat();
    ```
    
    
* Option 3: `then` style:

  ```chapel
  if foo > bar then doThis();
  else doThat();
  ```

## Functions

### Declarations

When breaking a function declaration across multiple lines (due to arguments, where clause(s), and/or other decorators), 
putting the open brace on its own line will improve readability

Good:

```chpl
proc foo(x) where reallyLongConditional1(x) && 
  reallyLongConditional2(x) && reallyLongConditional3(x) 
{
  ...
}
```

Bad
```chpl
proc foo(x) where reallyLongConditional1(x) && 
  reallyLongConditional2(x) && reallyLongConditional3(x) {
  ...
}
```


## Indexing

Brackets and parentheses are semantically equivalent for indexing in Chapel. Brackets are generally preferred over parentheses. 
Parentheses are acceptable in cases where indexing into a tuple or a custom type that defines a `this` method:

```chpl
var A: [1..10] int;
// OK
A[1];
// Not OK
A(1);

var t = (1,2,3);
// OK
t[2];
// Also OK
t(2);

record R {
  proc this(i) {
    return i + 1;
  }
}

var r = new R();
// OK
r[10];
// Also OK
r(10);
```

## Comments

### chpldoc comments

### non-chpldoc comments

For non-`chpldoc` comments, there a few styles, used in a variety of situations:

```chpl
//
// Clean and distinguishable comment usually used to label sections of code
//

foo();
```

```chpl
// Comment about the next line
foo();
```

This style is not used frequently, but is useful in some scenarios:
```chpl
foo(); // Comment about the current line
```n
