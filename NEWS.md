<!--
This file should document all pull requests and all user-visible changes.

When a pull request is completed, changes made should be added to a section at
the top of this file called "# UNRELEASED". All changes should be categorized
under "## MAJOR CHANGES", "## MINOR CHANGES", or "## BUG FIXES" following the
major.minor.patch structure of semantic versioning. When applicable, entries
should include direct links to the relevant pull requests.

Then, when a new release is made, "# UNRELEASED" should be replaced by a heading
with the new version number, such as "# CHANGES IN BMLePhoto VERSION 2.0.0."
This section will combine all of the release notes from all of the pull requests
merged in since the previous release.

Subsequent commits will then include a new "UNRELEASED" section in preparation
for the next release.
-->

# UNRELEASED

# CHANGES IN BMLePhoto VERSION 1.0.6

## BUG FIXES

- Avoid repeated ePhoto ODE integrations by caching complete results within each
  leaf conductance solve. A cached result is reused only when all four inputs
  match exactly: leaf temperature (`Tleaf`), incident PPFD (`ephoto_ppfd`),
  intercellular CO2 concentration (`Ci`), and experiment identifier (`exp_id`).
  There is no rounding or interpolation. Assimilation (`A`) is a cached output,
  not a cache input. Results are not shared across leaves or timesteps.
- Use the previous converged Ci and a local residual slope to guide subsequent
  bracket searches during leaf conductance coupling. Retain bracket and
  convergence checks, with a fallback to the original starting estimate when
  the previous-solution search fails.
- Restore complete cached biochemical outputs at the accepted Ci root while
  recomputing conductance for the current inputs. Biochemical equations and
  solver tolerances are unchanged.

# BMLePhoto VERSION 1.0.0

- This is the initial release of the package.
