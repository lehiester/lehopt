## "Textbook" simplex implementation

The "textbook" implementations are not modeled after any one book or course but
rather attempt to use what seems to be the most common options that are *fully
defined*; brief allusions to more sophisticated options without elaborating on
critical implementation details are common but here are not considered adequate
for inclusion in these implementations.

Detailed expositions of dual simplex are, of course, much less common; the
guiding principle for the textbook dual simplex implementation here was to have
the clearest parallels to textbook primal simplex as defined above. It is
mainly motivated by inclusion in a "textbok" branch-and-bound implementation
for which the ubiquity of dual simplex and warm-start is here considered
sufficiently well known.

The following implementation choices have been made:
- Standard form (equality constraints with nonnegative variables)
- No basis matrix factorization updating; recalculated every iteration
- No updating of primal values and dual slack; recalculated every iteration
  - Non-textbook critical details: What to do about error accumulation? What's
    the appropriate recalculation interval/rule to use?
  - Should not be performance bottleneck without factorization updating anyway
