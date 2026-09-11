# Current production Bias40 cone frame

The production analytic branch uses the authoritative local pipe tangent `t` at the evaluation arc position. It constructs a projected global-`+Z` normal

`n = normalize(+Z - (+Z·t)t)`

and binormal `b = normalize(t × n)`. With `cone-axis-bias-deg=40`, the biased cone axis is

`c0 = normalize(cos(40°) t + sin(40°) n)`.

The transverse basis is `e10 = normalize(n - (n·c0)c0)` and `e20 = normalize(c0 × e10)`. The field direction is the production formula

`Bhat = cos(30°) c0 + sin(30°)[cos(ψ)e10 + sin(ψ)e20]`, with `ψ(t)=ψ0+sense·2π·30t`.

`χ=0°` is exactly this legacy frame. A nonzero `cone-frame-azimuth-deg=χ` applies one Rodrigues rotation about the same local tangent `t` to **all three** vectors (`c0`, `e10`, `e20`) before the unchanged field formula is evaluated. It is therefore a spatial cross-sectional azimuth, not a phase offset.

The construction preserves orthonormality, handedness, the 40° axis-to-tangent bias, and the 30° cone half-angle to numerical precision. Existing jobs omit the option and retain the legacy frame.
