# Maintenance conventions: a gradual refactoring pilot

ctmaPower is the pilot. Other public ctma functions retain their interfaces and
implementations. Apply the following conventions when extending them; do not
rewrite fitting engines and statistical interpretation simultaneously.

1. Public functions orchestrate validation, fitting, calculations and presentation.
   Keep named arguments/backends at this boundary. Never silently ignore controls.
2. Pure numerical helpers accept explicit matrices/designs and return named
   results. They do not fit CT models, inspect global options or save files.
3. Extract backend/legacy fit data through one coordinate-aware adapter. Record
   matrix name, outcome row, predictor column, units and fixed/free status rather
   than guessing parameter order from labels or recycled vectors.
4. Separate effect calculation from numerical probability/sample-size inversion.
   Document the estimand, assumptions, reference probability convention, boundary
   behavior and missingness policy next to the implementation.
5. Keep file I/O explicit. Caches must identify their model, covariance/design and
   reference sample size. Do not treat an old filename as proof of compatibility.
6. Print concise stage messages. Retain method decisions and interpretation notes
   in returned summaries so they survive saving/reloading and can be audited later.
7. Test mathematical results and public behavior. Include negative effects,
   asymmetric matrices, dimensions above two, two-wave inputs, fixed parameters,
   missing metadata, singular generators, extreme effects and tied optima.
8. Verify optimization-free calculations against an independent fitting route
   where the models are equivalent. Do not broaden equivalence claims to
   constrained/non-Gaussian/hierarchical models without separate justification.
9. Keep experimental backend adaptation separate from shared statistical helpers.
   CoTiMA and CoTiMA2 share this numerical module; their public fitting-control
   adapters differ. Changes to shared helpers require tests in both packages.
10. Preserve old code through Git instead of backup R files inside R/. Build from
    clean source snapshots and avoid replacing unrelated working-tree changes.

## ctmaPower module map

* Public ctmaPower: resolves backend controls, obtains the invariant fit, invokes
  the analysis and assembles the established result object.
* .ctmaPowerDesign / Gaussian / RegressionModel: validates the design and decides
  whether direct unrestricted Gaussian covariance regression is justified.
* .ctmaPowerMatrices / Implied: extracts coordinate-labelled matrices, restores
  original time units and calculates CT-implied covariance moments.
* .ctmaPowerDirect / Lavaan: computes each full/reduced regression effect once
  per lag, independently of requested power levels.
* .ctmaPowerProbability / RequiredN: evaluates a documented probability convention
  and searches integer sample sizes without sequential scanning.
* .ctmaPowerOptima / Result: reports actual minimizing lags and structured output.

Do not infer a DT zero coefficient from a CT zero: indirect paths can invalidate
that inference. The default generated DT regression is unrestricted; explicit
DT restrictions belong in regressionModel. Gaussian metadata alone is also
insufficient: model structure must pass the unrestricted-regression check.
