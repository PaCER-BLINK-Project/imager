# Optional FRB scattering

The visibility injector can broaden each intrinsic boxcar pulse with a normalized,
causal exponential impulse response. The original 13-column catalogue is unchanged.
Append two or three optional fields to each event (whitespace-separated):

```
# obs_id id ra_deg dec_deg arrival_unix_s dm width_s peak_jy ref_mhz units_per_jy phase_sign tsub fsub [scattering_tau_s scattering_ref_mhz [scattering_index]]
1276619416 demo 92.0 -13.0 1592584200.25 100 0.080 10 170.36 1 -1 0 32 0.050 150 -4.0
```

- `scattering_tau_s`: exponential e-folding timescale in seconds, nonnegative.
  Zero (default) disables scattering.
- `scattering_ref_mhz`: finite positive reference frequency in MHz. Required
  if the scattering timescale is supplied in the catalogue. For C++ the
  field is `scattering_ref_hz` (default 150e6 Hz).
- `scattering_index`: exponent in `tau(nu) = tau_ref * (nu / nu_ref)^index`,
  finite, defaults to -4.0 if omitted. Positive indices and zero are supported.

The dispersion reference (`ref_mhz`) is independent of the scattering reference.
`arrival_unix_s` still refers to the dispersed leading edge at `ref_mhz`, not the
peak of the scattered pulse. `width_s` and `peak_jy` describe the **intrinsic**
boxcar, before scattering. The intrinsic fluence is `peak_jy * width_s` Jy s.

For time x since the dispersed leading edge, intrinsic amplitude A and width W,
the scattered flux density is

```
S(x) = 0                                            (x < 0)
     = A * (1 - exp(-x/tau))                         (0 <= x <= W)
     = A * (1 - exp(-W/tau)) * exp(-(x-W)/tau)        (x > W)
```

The implementation uses stable `expm1` expressions and analytic integration
within each time bin when `tsub=0` (recommended). The original midpoint sampling
is available when `tsub>0`. Frequency averaging retains the existing `fsub`
midpoint samples, with dispersion and scattering evaluated independently at
each subfrequency. The original unscattered behaviour is retained at tau=0.

The profile is truncated at `W + 25*tau` at each subfrequency, losing at most
`exp(-25)` (~1.39e-11) of the intrinsic fluence before temporal/frequency
quadrature and float conversion. No renormalisation is applied. The sparse
injection plan includes the tail even when the intrinsic pulse began in an
earlier visibility buffer. Computational cost increases with tail duration.
The implementation assumes temporal broadening only (an unresolved source),
not a scattering disk, scintillation, or a frequency-dependent intrinsic spectrum.

Use the same catalogue and absolute reference time across all independently
processed bands and buffers. In C++:

```cpp
// After constructing an existing blink_frb::Event:
event.scattering_tau_s = 0.050;
event.scattering_ref_hz = 150e6;
event.scattering_index = -4.0;
```

The existing `frb_injection_test` now checks the analytical profile against
independent quadrature, fluence conservation, frequency scaling, midpoint
sampling, catalogue parsing, split-buffer parity, and visibility injection
through the AstroIO adapter. Run it with the existing CTest configuration:

```
ctest --test-dir build -R '^frb_' --output-on-failure
```
