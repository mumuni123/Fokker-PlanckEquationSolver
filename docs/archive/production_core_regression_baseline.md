# Production Core Regression Baseline

This baseline freezes the accepted production reference summarized in
`18fs正式长跑问题与根因分析.md` section 15.1.  It is an acceptance gate for
subsequent algorithm changes, not a request to fit a different solver point by
point.

| time | core max abs(Ex)/E0 | core background density | beam front | beam max density |
|---:|---:|---:|---:|---:|
| 12 fs | 6.36 | 0.938-1.050 | 3.58 um | 1.54 nb0 |
| 15 fs | 15.39 | 0.871-1.149 | 4.47 um | 8.86 nb0 |
| 18 fs | 23.16 | 0.818-1.258 | 5.36 um | 8.03 nb0 |

At 18 fs the total-current range is `[-1.406, 0.914]e18 A/m2`.

For each checkpoint, retain the production snapshots for `Ex`, background
density, beam density, and total current.  Compare the core interval
`0.2-5.6 um` using normalized L2 errors, main peak and wave-front positions,
zero crossings or wavelength, beam-front and bunching-peak positions, and
total-current positive/negative peaks.

The first-level gate is:

- core field peak change <= 5 percent;
- field peak displacement <= 0.1 um;
- beam-front displacement <= 0.05 um;
- 18 fs total-current extrema change <= 10 percent;
- no more than 10 percent global suppression of core background-density
  oscillation; and
- no sudden multi-fold increase in the core FCT activation rate.

Boundary-only spikes are excluded from this core gate.  The gate must not be
satisfied by changing `J_bkg`, applying a Poisson/zero-mode projection, global
distribution scaling, global FCT strengthening, or replacing the core
high-order flux with first-order upwind transport.
