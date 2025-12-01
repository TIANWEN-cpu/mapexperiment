"""Lambert conformal conic projection utilities.

The script implements the two standard parallels version of the Lambert
conformal conic projection. It follows the steps shown in the supplied
course material: compute the basic helper values *m* and *t*, derive the
projection constants *n*, *F* (denoted as ``K`` in some texts), the radius of
curvature ``rho`` for any latitude, and finally project geographic
coordinates to the planar grid.

The implementation is self contained and exposes both a Python API and a
minimal command line interface for quick experiments.
"""
from __future__ import annotations

import argparse
import dataclasses
import math

from typing import Iterable, Tuple


@dataclasses.dataclass(frozen=True)
class Ellipsoid:
    """Representation of a reference ellipsoid.

    Attributes
    ----------
    a:
        Semi-major axis in meters.
    f:
        Flattening of the ellipsoid.
    """

    a: float
    f: float

    @property
    def eccentricity(self) -> float:
        """First eccentricity of the ellipsoid."""

        return math.sqrt(2 * self.f - self.f * self.f)


@dataclasses.dataclass(frozen=True)
class LambertParameters:
    """Projection constants for a Lambert conformal conic projection."""

    ellipsoid: Ellipsoid
    phi1: float
    phi2: float
    phi0: float
    lambda0: float
    n: float
    F: float
    rho0: float


# Default ellipsoid: WGS84
WGS84 = Ellipsoid(a=6_378_137.0, f=1 / 298.257_223_563)


def deg_to_rad(angle_deg: float) -> float:
    """Convert degrees to radians."""

    return math.radians(angle_deg)


def rad_to_deg(angle_rad: float) -> float:
    """Convert radians to degrees."""

    return math.degrees(angle_rad)


def _m(phi: float, eccentricity: float) -> float:
    """Helper function *m* from the textbook.

    m = cos(phi) / sqrt(1 - e^2 * sin^2(phi))
    """

    sin_phi = math.sin(phi)
    return math.cos(phi) / math.sqrt(1 - (eccentricity**2) * (sin_phi**2))


def _t(phi: float, eccentricity: float) -> float:
    """Helper function *t* from the textbook.

    t = tan(pi/4 - phi/2) * ((1 - e * sin(phi)) / (1 + e * sin(phi))) ** (e/2)
    """

    sin_phi = math.sin(phi)
    numerator = 1 - eccentricity * sin_phi
    denominator = 1 + eccentricity * sin_phi
    ratio = numerator / denominator
    # The textbook defines t with the inverse of the eccentricity factor compared
    # to the Mercator expression: divide by ( (1 - e sin φ) / (1 + e sin φ) )^(e/2).
    return math.tan(math.pi / 4 - phi / 2) * (ratio ** (-eccentricity / 2))


def compute_parameters(
    phi1_deg: float,
    phi2_deg: float,
    phi0_deg: float,
    lambda0_deg: float,
    ellipsoid: Ellipsoid = WGS84,
) -> LambertParameters:
    """Derive the Lambert projection constants.

    Parameters follow the notation from the course notes:
    - ``phi1`` and ``phi2``: standard parallels in degrees.
    - ``phi0``: latitude of origin in degrees.
    - ``lambda0``: central meridian in degrees.
    """

    e = ellipsoid.eccentricity
    phi1 = deg_to_rad(phi1_deg)
    phi2 = deg_to_rad(phi2_deg)
    phi0 = deg_to_rad(phi0_deg)
    lambda0 = deg_to_rad(lambda0_deg)

    m1, m2 = _m(phi1, e), _m(phi2, e)
    t1, t2 = _t(phi1, e), _t(phi2, e)

    n = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
    F = m1 / (n * (t1**n))

    t0 = _t(phi0, e)
    rho0 = ellipsoid.a * F * (t0**n)

    return LambertParameters(ellipsoid, phi1, phi2, phi0, lambda0, n, F, rho0)


def project(lon_lat_deg: Iterable[Tuple[float, float]], params: LambertParameters):
    """Project a sequence of ``(lon, lat)`` coordinate pairs.

    Returns a list of ``(x, y)`` tuples in meters.
    """

    e = params.ellipsoid.eccentricity
    results = []
    for lon_deg, lat_deg in lon_lat_deg:
        phi = deg_to_rad(lat_deg)
        lam = deg_to_rad(lon_deg)
        t = _t(phi, e)
        rho = params.ellipsoid.a * params.F * (t**params.n)
        theta = params.n * (lam - params.lambda0)
        x = rho * math.sin(theta)
        y = params.rho0 - rho * math.cos(theta)
        results.append((x, y))
    return results


def inverse_project(xy_pairs: Iterable[Tuple[float, float]], params: LambertParameters):
    """Invert Lambert conformal conic coordinates back to lon/lat."""

    e = params.ellipsoid.eccentricity
    results = []
    for x, y in xy_pairs:
        rho = math.hypot(x, params.rho0 - y)
        if params.n < 0:
            rho = -rho
        t = (rho / (params.ellipsoid.a * params.F)) ** (1 / params.n)

        # Iteratively solve for phi using Newton-Raphson
        phi = math.pi / 2 - 2 * math.atan(t)
        for _ in range(15):
            sin_phi = math.sin(phi)
            part = ((1 - e * sin_phi) / (1 + e * sin_phi)) ** (e / 2)
            phi_next = math.pi / 2 - 2 * math.atan(t * part)
            if abs(phi_next - phi) < 1e-12:
                phi = phi_next
                break
            phi = phi_next

        theta = math.atan2(x, params.rho0 - y)
        lam = params.lambda0 + theta / params.n
        results.append((rad_to_deg(lam), rad_to_deg(phi)))
    return results



def run_example():
    """Run an example mirroring the worked steps in the scanned notes."""

    # Parameters follow the table on page 107 (φ1=15°, φ2=50°, φ0=40°, λ0=105°E)
    params = compute_parameters(15.0, 50.0, 40.0, 105.0)

    sample_lon_lat = [(113.0, 40.0)]
    projected = project(sample_lon_lat, params)
    inverted = inverse_project(projected, params)

    print("Lambert conformal conic (two standard parallels)")
    print(f"Ellipsoid: a={params.ellipsoid.a:.3f} m, f={params.ellipsoid.f:.12f}")
    print(
        f"Parallels: φ1=15°, φ2=50°, origin φ0=40°, central meridian λ0=105°; n={params.n:.6f}, F={params.F:.6f}"
    )
    for (lon, lat), (x, y), (lon_inv, lat_inv) in zip(sample_lon_lat, projected, inverted):
        print(f"Input lon/lat : ({lon:.4f}°, {lat:.4f}°)")
        print(f"Projected     : x={x:,.3f} m, y={y:,.3f} m")
        print(f"Inverse check : ({lon_inv:.4f}°, {lat_inv:.4f}°)\n")



def _parse_pair(pair_str: str) -> Tuple[float, float]:
    try:
        lon_str, lat_str = pair_str.split(",")
        return float(lon_str), float(lat_str)
    except ValueError as exc:  # pragma: no cover - CLI convenience
        raise argparse.ArgumentTypeError(
            "Coordinate pairs must be specified as 'lon,lat'"
        ) from exc


def main():  # pragma: no cover - exercised via CLI
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parallels",
        nargs=2,
        type=float,
        metavar=("PHI1", "PHI2"),
        default=(15.0, 50.0),
        help="Standard parallels φ1 and φ2 in degrees (default: 15 50)",
    )
    parser.add_argument(
        "--origin-lat",
        type=float,
        default=40.0,
        help="Latitude of origin φ0 in degrees (default: 40)",
    )
    parser.add_argument(
        "--central-meridian",
        type=float,
        default=105.0,
        help="Central meridian λ0 in degrees (default: 105)",
    )
    parser.add_argument(
        "--coord",
        action="append",
        type=_parse_pair,
        metavar="LON,LAT",
        help="Longitude/latitude pair(s) to project (can be repeated)",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="Run the built-in example that mirrors the textbook values.",
    )

    args = parser.parse_args()

    if args.example:
        run_example()
        return

    coords = args.coord or [(113.0, 40.0)]
    phi1, phi2 = args.parallels
    params = compute_parameters(phi1, phi2, args.origin_lat, args.central_meridian)

    print("Projection constants")
    print(f"n = {params.n:.8f}")
    print(f"F = {params.F:.8f}")
    print(f"rho0 = {params.rho0:.3f} m")
    print()


if __name__ == "__main__":  # pragma: no cover
    main()
