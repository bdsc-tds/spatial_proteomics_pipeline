#!/usr/bin/env python3

import argparse
import colorsys
import json
import numpy as np
from functools import reduce
from pathlib import Path
from shapely.geometry import MultiPolygon, Polygon, mapping
from shapely.geometry.polygon import orient
from skimage import measure

from skimage import measure, filters
import numpy as np

def extract_contours_smooth(crop, r0, c0):
    """
    Smooth contour extraction that truly follows the mask.
    Blurs the binary mask slightly so the 0.5 isocontour is smooth,
    without ever moving outside the original mask boundary.
    """
    # Pad to avoid edge artifacts
    pad = 2
    padded = np.pad(crop.astype(float), pad)
    
    # Gaussian blur on the binary mask — the isocontour at 0.5
    # stays inside the mask but traces it smoothly
    blurred = filters.gaussian(padded, sigma=1.5)
    
    # find_contours on the blurred field at 0.5 = true mask boundary
    contours = measure.find_contours(blurred, level=0.5)
    
    # Adjust offset to account for padding
    adjusted = []
    for c in contours:
        shifted = [[pt[1] + c0 - pad, pt[0] + r0 - pad] for pt in c]
        if len(shifted) < 4:
            continue
        if shifted[0] != shifted[-1]:
            shifted.append(shifted[0])
        adjusted.append([[float(p[0]), float(p[1])] for p in shifted])
    
    return adjusted

def rings_to_geometry(rings):
    """
    GeoJSON geometry of one label from its closed contour rings.
    Rings nested inside another ring are holes (even-odd rule), so e.g. a cytoplasm ring
    becomes one polygon with the nucleus as a hole instead of two filled polygons.
    """
    # the 0.5 isolines never cross, so XOR-ing the ring polygons yields the even-odd fill
    polygons = [Polygon(r) for r in rings]
    polygons = [p if p.is_valid else p.buffer(0) for p in polygons]
    geometry = reduce(lambda a, b: a.symmetric_difference(b), polygons)

    # counter-clockwise exteriors and clockwise holes, as RFC 7946 recommends
    if isinstance(geometry, Polygon):
        geometry = orient(geometry)
    elif isinstance(geometry, MultiPolygon):
        geometry = MultiPolygon([orient(p) for p in geometry.geoms])
    else:
        return None
    return None if geometry.is_empty else mapping(geometry)

def label_to_hex(label: int) -> str:
    golden = 0.618033988749895
    hue = (label * golden) % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.95)
    return "#{:02X}{:02X}{:02X}".format(int(r * 255), int(g * 255), int(b * 255))

def contour_to_geojson_coords(contour, row_offset, col_offset):
    coords = [[float(pt[1] + col_offset), float(pt[0] + row_offset)] for pt in contour]
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    return coords

def mask_to_geojson(npy_path, output_path=None, min_pixels=5):
    mask = np.load(Path(npy_path))
    assert mask.ndim == 2, f"Expected 2D mask, got {mask.shape}"
    print(f"Mask shape: {mask.shape}, unique labels: {len(np.unique(mask))}")

    props = measure.regionprops(mask)
    features = []
    skipped = 0

    for region in props:
        if region.area < min_pixels:
            skipped += 1
            continue

        r0, c0, r1, c1 = region.bbox
        crop = (mask[r0:r1, c0:c1] == region.label).astype(np.uint8)

        # Pad if too small
        ph = max(0, 2 - crop.shape[0])
        pw = max(0, 2 - crop.shape[1])
        if ph > 0 or pw > 0:
            crop = np.pad(crop, ((0, ph), (0, pw)))

        contours = measure.find_contours(crop, level=0.5)
        if not contours:
            skipped += 1
            continue

        contours_sorted = sorted(contours, key=len, reverse=True)
        rings = extract_contours_smooth(crop, r0, c0)
        if not rings:
            skipped += 1
            continue

        # polygon(s) with holes from the rings
        geometry = rings_to_geometry(rings)
        if geometry is None:
            skipped += 1
            continue

        features.append({
            "type": "Feature",
            "geometry": geometry,
            "properties": {
                "objectType": "cell",
                "color_hex": label_to_hex(int(region.label)),
                "label": int(region.label),
            },
        })

    print(f"Features: {len(features)}, skipped: {skipped}")
    geojson = {"type": "FeatureCollection", "features": features}

    if output_path:
        with open(output_path, "w") as f:
            json.dump(geojson, f)
        print(f"Saved → {output_path}")

    return geojson


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert a .npy segmentation mask to GeoJSON."
    )
    parser.add_argument("npy",    help="Input .npy mask file")
    parser.add_argument("output", help="Output .geojson file")
    parser.add_argument("--min-pixels",  type=int,   default=5,   metavar="N",
                        help="Skip labels with fewer than N pixels (default 5)")
    args = parser.parse_args()

    mask_to_geojson(
        args.npy, args.output,
        min_pixels=args.min_pixels,
    )


if __name__ == "__main__":
    main()