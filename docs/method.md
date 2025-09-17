# SOC scaling

The workflow scales the GFN1-xTB atomic SOC constants in DFTB+ to fit
selected periodic spin–orbit splittings, with one scale per configured
element and shell. The constants come from
[Jha and Heine, Table 1](https://doi.org/10.1021/acs.jctc.2c00376).

## Samples

DFTB+ and GPAW use the same primitive cell for each sample. For the
diamond and zincblende builders, an isotropic scalar GFN1-xTB relaxation
provides the lattice and the crystal symmetry fixes the fractional basis.
DFTB+ 24.1 uses the `ConjugateGradient` driver for constrained lattice
relaxation.
A material can instead supply a `structure_path` for an already prepared
parent cell, resolved relative to the config. Strained samples scale the
periodic axes of the parent at fixed fractional coordinates.

## Observables

Each observable targets a scalar manifold of size 1, 2 or 3 at a special
k-point, which must split into spinor multiplets of sizes 1+1, 2+2 or 2+4.
Splittings are differences between the mean energies of the two
multiplets. Degeneracy uses the configured energy-span tolerance.
Conduction selection starts at the first degenerate manifold above the
scalar gap.

Each observable also records the scalar multiplicities from the gap to the
target. For example, `1,2` identifies a conduction doublet preceded by a
singlet. Only equal patterns enter the fit. A mismatched row receives the
`manifold_mismatch` flag.

For matching patterns, let $\Delta_{\mathrm{xTB}}$ and $\Delta_{\mathrm{ref}}$
be the gaps to the closest neighbouring scalar manifold in the two spectra.
The row is flagged `unisolated_manifold` when
$|y_{\mathrm{ref}}|\,|1/\Delta_{\mathrm{xTB}}-1/\Delta_{\mathrm{ref}}|>0.1$.
The test uses the L-1 xTB scalar spectrum and the all-tight GPAW
scalar spectrum. Orbital character remains untested. The signed SOC-induced shift of the
mesh band gap, ΔE_g, is a diagnostic.

## Reference and σ

GPAW 25.7.0 converges a non-spin-polarized scalar PW/PBE density on a
Monkhorst–Pack mesh. The fixed-density step computes scalar eigenstates at
the same special k-points used in DFTB+ with symmetry disabled and with
every retained band converged. GPAW diagonalizes the SOC Hamiltonian in
that scalar band basis after SCF. DFTB+ runs SOC self-consistently in SCC.

The config sets base and tight values for the cutoff, mesh and scalar band
count. Three one-knob runs tighten each of them individually. The target
comes from the all-tight run. σ is the largest absolute difference of an
individual one-knob run from the base.

The `reference` config field also accepts a path to a CSV, resolved
relative to the config. It must contain one row per applicable sample and
configured observable. Its columns are
`sample,observable,value_ev,sigma_ev,delta_ecut_ev,delta_kmesh_ev,delta_nbands_ev,pattern,neighbour_gap_ev,source`.
`neighbour_gap_ev` is the positive gap to the closest scalar manifold.
The computed CSV takes it from the all-tight run.
Use comma-separated multiplicities in the `pattern` field and quote it
when it contains a comma.

## Loss, prior and identifiability

The data residual is
$r_i = (y_{\mathrm{xTB},i}-y_{\mathrm{ref},i})/[\eta\max(|y_{\mathrm{ref},i}|,\sigma_i)]$,
where $\eta$ is the configured `relative_error`. The objective includes one
Gaussian prior residual $(s_p-1)/\tau$ per active scale, with
$s_p \in [0.5,1.5]$. Reduced χ² is $\sum_i r_i^2/(n-p)$, where $p$ is the
number of active scales and $n$ the number of observations. It is reported
when $n>p$.

At scale 1 the finite-difference Jacobian estimates the change across the
full bound width. A scale is identifiable if its change exceeds 3σ for at
least one eligible observable. It tests each scale individually.
Unidentifiable scales stay at 1. The optimizer uses SciPy trust-region
reflective least squares. Scale standard errors are the square roots of
the diagonal of $(J^\mathsf{T}J)^{-1}$, with prior rows in $J$.

## L-1 through L2

- **L-1** uses the original constants at scale 1.
- **L0** fits one scale per element and shell across the training samples.
  A valence-only L0 control is also evaluated on the eligible conduction rows.
- **L0 common** fits a constant scale per element and shell in each training
  fold using the same observables as L1.
- **L1** fits one scale per sample and shell using the observables that
  pass the checks for every sample of that material.
- **L2** fits ridge regression to the natural log of identifiable L1
  scales, separately for each material and species. The descriptor is
  the ASE neighbour-list mean of the nearest-neighbour distances of the
  fitted species. Descriptor statistics are fit on the training fold.
  The weight is the inverse squared log-scale standard error, estimated
  as $\operatorname{SE}(s)/s$. The intercept is unpenalized. Predicted
  scales are bounded to the fitting interval.

Each leave-one-sample-out fold fits L0, L0 common and L2 using the remaining samples.
DFTB+ is evaluated at each predicted held-out scale and compared to the
held-out GPAW splittings. The MAE for all four levels uses the same set
of eligible rows. The trend table fits an ordinary least-squares slope to
each splitting against distance across all samples, using the held-out
predictions for the three fitted levels and the same observables as L1.
Each fitted forward evaluation records the scales and predictions in
`evaluations.csv` alongside its DFTB+ inputs and outputs. The scalar
DFTB+ spectra are computed in `baseline` and reused throughout. Each fit
evaluation reconverges SOC from the corresponding L-1 SCC charges.
