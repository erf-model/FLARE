# Broken-wire state at rest

## Purpose

`broken_wire.py` finds the state at rest, in still air, of a conductor section: a row of spans
between two dead ends that hang from insulator strings at the towers between (or are clamped to
them). The section may be whole, or one span may have let go of one of its supports, or one
string may have left its cross-arm. The tool gives the tension left in the spans that stay, the
angle each string leans to, and the pull each tower and dead end then carries.

It uses the balance of forces only: every span and string is a heavy elastic cable, and at the
lower end of every string the pulls of the cables meeting there add up to zero. It shares no
code with ERF or MoorDyn. The conductor model's unit tests `BrokenWire.*` compare MoorDyn-C at
rest after a release (`conductors.failure_break`, `conductors.tower_failure`) with its numbers.
The equations are in "State at rest after a release" in `Docs/sphinx_doc/theory/Conductors.rst`.

It is Python with the standard library only and runs no model.

## Files

- `broken_wire.py`: the tool.
- `tests/test_broken_wire.py`: its tests (17, in the classes Cable, Balance, Limits, Robust, Reference, Refused).
- `BrokenWireTests.cmake`: the CTest registrations, included by `Tests/CTestList.cmake`.

## Data and acknowledgements

No data set is used or redistributed, and no published formula is taken over: the cable's
equations are derived in the documentation from the balance of forces on a piece of cable and
Hooke's law. The level span is checked against the elastic catenary the model's documentation
already uses (its "Verification" section). The tool does not reproduce a design standard's
broken-wire load; none was available to check against.

## The plan

A JSON object. Lengths in m, forces in N. The section lies in one vertical plane: `x` along the
line, `z` up.

| Key | Default | Range | Effect |
|---|---|---|---|
| `end_a`, `end_b` | required | `[x, z]`, finite | the dead ends |
| `towers` | none | list of `[x, z]`, `x` increasing from `end_a` to `end_b` | the cross-arms the strings hang from, or the clamps |
| `lengths` | required | one per span, > 0 | unstretched length of each span |
| `diameter` | required | (0, 1] | conductor diameter |
| `mass_per_length` | required | > 0, heavier than the air displaced | conductor mass (kg/m) |
| `axial_stiffness` | required | > 0 | conductor EA |
| `insulator_length` | 0 | >= 0 | string length; 0: the conductor is clamped at the towers |
| `insulator_mass` | required with strings | > 0 | mass of one string (kg) |
| `insulator_diameter` | 0.254 | > 0 | disc diameter, for the string's buoyancy |
| `air_density` | 1.2 | >= 0 | air density (kg/m^3), for buoyancy |
| `gravity` | 9.81 | > 0 | m/s^2 |
| `ground` | 0 | below every attachment | height of the flat ground a released span falls to |
| `segments` | 0 | integer in [0, 100000] | 0: continuous cables; n: MoorDyn's lumped chain of n links per span (2 per string) |
| `released` | none | `{"span": k, "end": "a" or "b"}` or `{"string": j}` | what has let go: span `k` (from 1) at its end towards `end_a` or `end_b`, or the string at tower `j` |

A key outside its range, or an unknown key, stops the tool with a message that names it (exit
code 2). Exit code 3: no state at rest was found. Exit code 4: a conductor still held at both
ends, or one that has lost its string, would lie on the ground, which the tool does not describe.

With `segments` a released span that reaches the ground loads its support in whole nodes. For a
height at which one of its nodes just touches the ground the chain has no state at rest and the
tool stops with exit code 3; the continuous cable (`segments` = 0) has one for every height.

## How to run

```bash
python3 broken_wire.py plan.json
```

```bash
python3 broken_wire.py plan.json --json
```

The first form prints a table; the second the same results as JSON.

## Expected results

For spans of 300, 500 and 100 m of 795 kcmil ACSR on 2.5 m strings of 60 kg, every attachment
60 m up, with span 1 released at tower 1:

```json
{"end_a": [100, 60], "end_b": [1000, 60], "towers": [[400, 60], [900, 60]],
 "lengths": [301.5, 506.944, 100.056], "diameter": 0.0281, "mass_per_length": 1.628,
 "axial_stiffness": 3.0e7, "insulator_length": 2.5, "insulator_mass": 60.0,
 "released": {"span": 1, "end": "b"}}
```

```
conductor weight 15.963379 N/m; continuous cables
span  attached  H (N)        T_a (N)      T_b (N)      lowest point above ground (m)
1     end a     hangs 59.985 m, 241.515 m on the ground; tension at the held end 957.563 N
2     both      11779.337    12467.976    12441.982    15.9717
3     both      11802.359    11813.070    11852.981    56.8276
tower angle (deg)  T_top (N)    pull on the tower: along the line, vertical (N)
1     69.6071      12672.541    11779.337    -4673.385
2     0.2750       5096.358     23.022       -5096.306
```

The whole section carries 13.57 kN in span 2 and pulls on tower 1 with 0.11 kN along the line.
After the release span 2 keeps 11.78 kN, 87 % of its tension, the string at tower 1 leans 69.6
degrees into it, and tower 1 carries that 11.78 kN along the line.

## CTest

`BrokenWire_Cable`, `BrokenWire_Balance`, `BrokenWire_Limits`, `BrokenWire_Robust`,
`BrokenWire_Reference` and `BrokenWire_Refused` (label `unit`), registered when `ERF_ENABLE_STRUCTURES` is ON and Python 3
is found. They run the tool only:

```bash
python3 tests/test_broken_wire.py
```
