# SOC fitting

This repository is a research prototype of a protocol that fits
shell-resolved spin–orbit coupling (SOC) scales for GFN1-xTB in DFTB+ to
periodic GPAW references. The open-source GPAW and DFTB+ stack keeps the
protocol reproducible; it is not intended for production parameterization.
The fitting examples cover diamond Si at five strains from scalar cell
relaxation through held-out evaluation, and an L-point fit for GaAs and
GaP. Their reference values are provisional.

See the [method](docs/method.md) and the result tables for
[Si strain](examples/results/si_strain_demo/report.md) and
[GaAs/GaP](examples/results/zb_lpoint_boundary/report.md).

## Run locally

The example environment uses Python 3.13.7,
GPAW 25.7.0, and DFTB+ 24.1.

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -c environment/python-constraints.txt \
    pip setuptools==80.9.0 wheel -r environment/python-wheels-osx-arm64.txt
.venv/bin/python -m pip install --no-build-isolation \
    -c environment/python-constraints.txt -e '.[dev,dft]'
conda create --prefix ./.dftb-24.1 --file environment/dftb-osx-arm64.txt
```

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
DFTBPLUS=./.dftb-24.1/bin/dftb+
for stage in relax reference baseline fit learn report; do
    .venv/bin/python -m ml_soc_xtb.periodic_workflow \
        --config configs/si_strain_demo.py --dftbplus-executable "$DFTBPLUS" \
        --jobs 5 "$stage" || break
done
```

The run writes `report.md`, `scales.png`, and the CSVs under
`examples/generated/si_strain_demo/`.
For the GaAs/GaP check, use [its config](configs/zb_lpoint_boundary.py) and run
`relax reference baseline fit report`.

Use `--output-root` to choose an output directory. Each directory stores the
original config bytes. Use a fresh directory when changing a config.

## Tests

```sh
.venv/bin/python -m pytest -q -m 'not dft_integration'
.venv/bin/python -m pytest -q -m dft_integration --dftbplus-executable="$DFTBPLUS"
```
