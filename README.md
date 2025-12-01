# Map projection experiment

This repository contains a small Python implementation of the Lambert conformal conic projection (two standard parallels) used in the map cartography exercise shown in the attached scans.

## Features
- Computes the projection constants `n`, `F` (often called `K` in textbooks), and `ρ0` from two standard parallels, a latitude of origin, and a central meridian.
- Projects geographic coordinates (longitude/latitude) to planar `(x, y)` values and performs the inverse transformation.
- Exposes both a Python API (`compute_parameters`, `project`, `inverse_project`) and a command line interface for quick calculations.

## Quick start
Run the built-in example that mirrors the values in the experiment (φ1=15°, φ2=50°, φ0=40°, λ0=105°) using the default WGS‑84 ellipsoid:

```bash
python main.py --example
```

Project custom coordinates with your own parameters:

```bash
python main.py \
  --parallels 20 45 \
  --origin-lat 25 \
  --central-meridian 110 \
  --coord 113,40 --coord 115,35
```

You can also import the module to reuse the functions programmatically:

```python
from main import compute_parameters, project
params = compute_parameters(20, 45, 25, 110)
xy = project([(113, 40)], params)
```
