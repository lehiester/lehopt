# TODO - immediate next steps for "textbook" MILP solver

- Add more slack variables for feasibility tolerance?
  - How to use the slack only "when needed"?
- gen-ip054 and gen-ip02 solutions are wrong! Their optima are not found
  - This must indicate the solver is deeming them to be within an
    infeasible node, right? Is there another explanation?
  - Can check if the optimal point is within each node and find which one it's deeming infeasible
- mad.mps appears to be cycling!
- B&B: include MIP gap in output
- B&B: best-first with memory-saving mode
- B&B node dual simplex: ability to pass in cutoff value
- Dual simplex: add degeneracy handling
- License notice for CSparse_modified



# TODO - Next in queue

- Invariant checking for SparseVector
  - Add checks for new invariant 01 (positive length!)
- Unit tests for SparseVector
  - Basic spot checks
  - Unit tests for panic conditions
- Unit tests for DenseVector
- Unit tests for CSCMatrix
  - (!) new_from_columns failing because of SparseVector::new bug!



# TODO

- Hunt down all uses of `.expect(fmt!(...))` and replace with
  `.unwrap_or_else(|| panic!(...))` to only call `fmt!` when panicking

- Rename and/or relocate "textbook" implementations of simplex and B&B

- For "non-textbook" stuff: status output on separate thread

- Replace all uses of `.try_into()` for `usize` values with new unchecked
  trait method which is a no-op when "converting" to `usize`
  - Make sure it's always checked against `Int::MAX` first!

- `Int::CHECK_CARDINALITY`: search for all uses of Int::MAX and wrap in an `if`
  - `I::MAX`
  - `CSInt::MAX`

- Review invariants for "fitting in I::MAX" and ensure that *length* fits, not
  the highest index! Highest index needs to be *less* than I::MAX so that I can
  also be used to represent the length! (Also revise verbiage)

- Unit tests specifically for limited-size-`Int` containers (smaller than
  `usize`), to ensure no mangling in general and proper application and error
  handling of size limits

- `Scalar`: require an "infinity" value
  - Add to trait definition
  - For `Rational`: zero denominator, positive numerator
      - Negative numerator for negative infinity
  - Replace use of `Option<S: Scalar>` in `MILPModel` with just a `Scalar` and
    the appropriate checking for if it `.is_finite()`, or something similar

- Result invalidation (e.g. do not allow querying solution by variable if
  variables may have been re-indexed)
  - Add an immutable reference to the `Result` struct? Then compiler would
    prohibit mutating the model until done using the result
  - Mutation nonce? Track the nonce that each variable was added?
  - Give each variable/constraint a persistent ID independent of its numeric
    index in the model data?
  - Probably a fair error: "attempted to read (look up, ...) deleted variable"
  - Differentiate between additive and destructive changes?



# Resolutions

- Design for petaflops-scale systems out of scope
  - Cardinality of sets indexed by 64-bit integers is not checked to ensure it
    fits within 64 bits, to reduce overhead
