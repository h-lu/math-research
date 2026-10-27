# Finite verification

Run `python3 repro/verify.py` from the manuscript directory, or `make check`.
Python 3 with the standard library is sufficient. No repository history,
external dataset, optional numerical package, or network connection is needed.

The program checks the mechanical counts, binary-prefix identities, source
and pullback degree inequalities, and the replicated Qin–Sun envelope.
Its default range is n = 2,...,20000, explicit monomial enumeration through
256, and target replication degrees 0,...,20000. It also checks the saved
exact fields in `data.json`; approximate display constants do not decide
pass or fail. `comparison.csv` is the finite comparison table.

`python3 repro/verify.py --table` prints the table. The default operation
writes no files. Use `--write-data` only to regenerate the saved dataset.
Run without `-O` or `-OO`, which the script explicitly rejects.

These are finite algebraic and counting checks. They do not independently
prove period independence, simple-zero persistence, simultaneous realization,
or the existence of the asserted limit cycles. Those arguments are in the
manuscript. No cubic-logarithmic lower bound is claimed.
