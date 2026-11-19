# sbml-io

Read and write [SBML](https://sbml.org) as plain Rust data. Pure Rust, so it
also works in WebAssembly.

> **Status: implemented, not yet published.** Lives in the PathwayDesigner
> workspace for now; it will move to its own repository when published.

## Scope

- **SBML Level 3**, Versions 1 and 2, all of core:
  - function definitions;
  - unit definitions;
  - compartments, species and parameters;
  - initial assignments;
  - rules: assignment, rate and algebraic;
  - constraints;
  - reactions with species references, modifiers and kinetic laws with local
    parameters;
  - events with trigger, delay, priority and event assignments.
- **SBML Level 2**, Versions 1 to 5, including what only Level 2 has:
  - compartment types and species types;
  - `charge`;
  - stoichiometry given as MathML;
  - a kinetic law's `listOfParameters`;
  - a unit's `offset`.
- **The Layout package** (Level 3 Layout Version 1), read and written in full.
  In Level 2 files, layout lives in an annotation. `sbml-io` reads it from there
  too, into the same `Layout` data.
- **The Render package** (Level 3 Render Version 1), read and written in full:
  render information (global and per layout), colours, gradients, line
  endings, styles, and every drawing primitive. In Level 2 files Render lives in
  annotations, and is read from there and written back there, into the same
  `render` data. Polygons written by older tools as `listOfCurveSegments` are
  read as points. Values that may be relative (`50%`, `5 + 50%`) are kept as
  written, with `RelAbs::parse` and `RelAbs::resolve` to read them.
- **MathML**, the subset SBML allows, as a tree.
- **Everything else is kept.** Notes and annotations, other packages (fbc,
  comp, …), and any attribute or element `sbml-io` doesn't understand all
  survive a read and a write unchanged.

Not in scope, at least at first:
- SBML Level 1;
- converting between levels and versions;
- checking units;
- the full set of libSBML's consistency rules;
- an infix formula syntax. That belongs to the application. MathML is offered as
  a tree to convert from and to.

## Design

**Not libSBML's object model.** There is no class hierarchy, no getters and
setters, no `ListOf` objects, no plugins and no casting. A document is a tree
of ordinary structs with public fields:
- `Vec` for lists;
- `Option` for attributes that may be absent;
- enums where SBML has alternatives: rules, MathML nodes, curve segments.

Building a model is `model.species.push(Species { … })`; reading one is
`for s in &model.species`.

```rust
pub fn read(xml: &str) -> Result<Read, Vec<Diagnostic>>; // errors carry line and column
pub fn write(doc: &Document) -> String;
pub fn validate(doc: &Document) -> Vec<Issue>;            // ids unique, references resolve

pub struct Read { pub document: Document, pub warnings: Vec<Diagnostic> }

pub struct Document {
    pub level: u32,
    pub version: u32,
    /// Namespace declarations on `<sbml>`, in order. Kept so that preserved
    /// XML from other packages stays valid.
    pub namespaces: Vec<(Option<String>, String)>,
    pub model: Option<Model>,
    pub meta: Meta,   // <sbml> may carry notes, annotations and package attributes
}

pub struct Model {
    pub id: Option<String>, pub name: Option<String>,
    // substanceUnits, timeUnits, … conversionFactor: Option<String>
    pub function_definitions: Vec<FunctionDefinition>,
    pub unit_definitions: Vec<UnitDefinition>,
    pub compartment_types: Vec<CompartmentType>,   // Level 2 only
    pub species_types: Vec<SpeciesType>,           // Level 2 only
    pub compartments: Vec<Compartment>,
    pub species: Vec<Species>,
    pub parameters: Vec<Parameter>,
    pub initial_assignments: Vec<InitialAssignment>,
    pub rules: Vec<Rule>,
    pub constraints: Vec<Constraint>,
    pub reactions: Vec<Reaction>,
    pub events: Vec<Event>,
    pub layouts: Vec<Layout>,                      // the Layout package
    pub layouts_meta: Meta,                        // listOfLayouts itself
    pub global_render: Option<render::RenderList>, // the Render package
    pub meta: Meta,
}

/// What any element may carry.
pub struct Meta {
    pub metaid: Option<String>,
    pub sbo_term: Option<String>,
    pub notes: Option<String>,       // the content of <notes>, raw XML
    pub annotation: Option<String>,  // the content of <annotation>, raw XML
    pub extra: Extra,
}

/// What `sbml-io` doesn't understand on an element, kept for writing back.
pub struct Extra {
    pub attributes: Vec<Attribute>,  // qualified name, namespace, value
    pub elements: Vec<String>,       // raw XML, in document order
}

pub enum Rule {
    Assignment { variable: String, math: Option<Math>, meta: Meta },
    Rate { variable: String, math: Option<Math>, meta: Meta },
    Algebraic { math: Option<Math>, meta: Meta },
}

pub enum Math {
    Number { value: Number, units: Option<String> },
    Ident(String), Symbol(Symbol), Constant(Constant),
    Apply(Op, Vec<Math>), Call(String, Vec<Math>), CallSymbol(Symbol, Vec<Math>),
    Root { degree: Option<Box<Math>>, arg: Box<Math> },
    Log { base: Option<Box<Math>>, arg: Box<Math> },
    Piecewise { pieces: Vec<(Math, Math)>, otherwise: Option<Box<Math>> },
    Lambda { args: Vec<String>, body: Box<Math> },
    Semantics { math: Box<Math>, annotations: Vec<String> },
    Raw(String),   // MathML outside the subset, kept as written
}

pub enum Segment {
    Line { start: Point, end: Point },
    Bezier { start: Point, end: Point, base1: Point, base2: Point },
}
```

**Rules that follow from this design:**
- **Attributes are stored as written.** Absent means `None`, not the default.
  Level 2 has defaults where Level 3 makes attributes required, and keeping what
  was written means a document read and written back doesn't gain attributes it
  didn't have. Applications that need the effective value apply the defaults of
  the document's level.
- **References between elements are id strings.** `validate` checks that they
  resolve.
- **A document is written at its own level and version.** The writer leaves out
  what is `None`, and writes kept XML back into the element it was found in.
- **Values are typed.** Numbers are `f64`, booleans `bool`. A number is
  written back the way SBML writers usually write it (`1`, `0.5`, `1e-20`,
  `INF`), not always as it was typed: `1.0` comes back as `1`.

**Where the implementation settled details the design left open:**
- Every `math` is an `Option<Math>`, because Level 3 Version 2 makes math
  optional.
- `Document` has a `meta` rather than only an `extra`: `<sbml>` itself may
  carry notes, annotations and package attributes such as `fbc:required`.
- Kept child elements are written after the children `sbml-io` knows, not at
  their exact original position. libSBML reads them either way.
- What can't be kept is reported as a warning: attributes on `listOf…`
  elements, and unknown elements inside them. Points, dimensions, bounding
  boxes and curves don't keep unknown XML.
- A raw piece of XML declares any namespace it inherited from an element
  below `<sbml>`, so that it stays valid wherever it is written.

## Testing

- Unit tests: `cargo test -p sbml-io`. The Layout and Render files in
  `tests/data` are written by python-libsbml
  (`tests/scripts/make_layout_fixtures.py`, `make_render_fixtures.py`); the
  Render ones use every part of the package, in Level 3 and in Level 2.
- Against libSBML: `examples/roundtrip.rs` reads and writes back every file in
  a folder, and `tests/scripts/libsbml_compare.py` has python-libsbml compare
  each rewrite with its original. CI runs this on every Level 2 and Level 3
  model of the [SBML Test Suite](https://github.com/sbmlteam/sbml-test-suite)
  (9,350 files in release 3.5.0): libSBML finds every rewrite identical to its
  original, including the parts from other packages that `sbml-io` keeps as
  raw XML. CI also runs it on the Layout and Render files in `tests/data`.

## Licence

MIT OR Apache-2.0.
