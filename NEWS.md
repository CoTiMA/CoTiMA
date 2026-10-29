# CoTiMA 1.0.4

* `ctmaTables()` now accepts `fontSize = 10` (points) for titles, headers,
  cells and notes in Word, XML, TeX and PDF exports.

## Regular Stan table compatibility (2026-10-29)

* Reads joint population DRIFT arrays from regular optimized and sampled Stan fits.
* Reconstructs missing optimized-fit caches from retained raw draws through the
  fitted Stan program, with Kalman filtering disabled and no model refitting.
* Correctly labels native ML, MAP and Bayesian fits from their direct args fields.
* Retains native-summary precision for four-decimal tables, avoids partial cache
  name matching, and preserves a wrapper's CT representation during DT refresh.
* Validated against regular ctsem 3.11.1 and an existing 2000-draw CoTiMA benchmark.



## Table presentation options (2026-10-29)

* layout selects portrait (long, default) or landscape (wide) pages consistently
  across Word, XML, TeX and PDF exports. Oversized tables continue in column panels.
* sortBy groups ctmaInit estimates by study (default) or by parameter.
* rowLayout keeps one parameter per row (default), combines DRIFT entries, or
  combines all selected parameters with grouped parameter/statistic headers.
* Numeric long-format values and safe staged output publication are preserved.



## Saved-fit table summaries (2026-10-29)

* ctmaTables defaults to four display decimals.
* Missing statistics and DT uncertainty can be regenerated from saved native fits
  with summary(timeinterval=lag), with one consent for all primary studies.
* recomputeSummary="yes" permits batch use; no models are refitted and raw point
  estimates are never treated as uncertainty draws.



## Table selectors and primary-study summaries (2026-10-28)

* Accepts bare, case-insensitive matrix families such as parameters=c(drift, diffusion).
  Quoted selectors and character-vector variables remain supported.
* Presents every primary study from a ctmaInit aggregate summary with original study
  labels, retained uncertainty and explicit aggregate fit-statistic notes; no refitting.
* Preserves unavailable statistics as unavailable and aligns study/parameter text
  columns in Word and TeX exports.



## APA tables (2026-10-28)

* Adds ctmaTables() with checked summary agreement, named/matrix parameter selection,
  source-aware statistics, interactive selection and APA-style Word/XML/TeX/PDF exports.
* Adds joint-draw DT effects at selected fitted-time lags and coordinate-aligned model comparisons.
* Exports are staged; existing outputs require explicit overwrite. Unsupported summaries
  and unavailable statistics fail explicitly rather than fabricating transformed values.
* Documents optional Word/PDF dependencies and the limits of legacy summary provenance.



## Power analysis refactoring (2026-10-28)

* Separates validation, invariant fitting/loading, matrix extraction, implied moments,
  regression effects, power probabilities, integer sample-size search and result assembly.
* Verifies unrestricted Gaussian regression eligibility for direct covariance calculations.
  Supports explicit restricted Gaussian regression syntax through a lavaan fallback, with
  console and stored summary notes identifying the calculation route.
* Corrects negative-effect p-values, two-wave and higher-dimensional calculations, lag
  bounds, reduced hypotheses, root boundaries and tied optimal-lag selection.
* Fast default power includes both rejection tails. useMBESS=TRUE retains the legacy
  probability convention. Missing post-hoc powers stay NA instead of becoming zero.
* Restores original time units from invariant-fit scaling metadata. Reduced-fit caches
  now validate covariance/model metadata and are reused across requested power levels.



## Plotting corrections (2026-10-28)

* Correct forest pooled confidence intervals and funnel confidence limits.
* Correct time-range indexing, axis positions, power-only plots and fractional-lag selection.
* Handle multiple fits, legacy example objects and per-fit moderator/lag mappings.
* Export uniquely named PNG files in the requested directory, including headless sessions.

