#!/usr/bin/env python3
import argparse
import numpy as np
import tifffile as tiff


def to_uint8(a: np.ndarray) -> np.ndarray:
    if a.dtype == np.uint8:
        return a
    if np.issubdtype(a.dtype, np.integer):
        return np.clip(a, 0, 255).astype(np.uint8)
    a = np.nan_to_num(a)
    mx = float(np.max(a)) if a.size else 0.0
    if mx <= 1.0:
        a = a * 255.0
    return np.clip(a, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser(description="Convert 3x8-bit grayscale TIFF/OME-TIFF to interleaved RGB TIFF.")
    ap.add_argument("input_tif")
    ap.add_argument("output_tif")
    ap.add_argument("--compression", default="deflate",
                    help="Compression: deflate, lzw, zstd, none. Default: deflate")
    ap.add_argument("--bigtiff", action="store_true", help="Force BigTIFF output.")
    args = ap.parse_args()

    compression = None if args.compression.lower() in ("none", "no", "false", "0") else args.compression

    img = tiff.imread(args.input_tif)

    # Normalize to (H, W, 3)
    if img.ndim != 3:
        raise ValueError(f"Expected 3D array, got shape={img.shape}")
    if img.shape[-1] == 3:
        rgb = img
    elif img.shape[0] == 3:
        rgb = np.transpose(img, (1, 2, 0))
    else:
        raise ValueError(f"Expected (H,W,3) or (3,H,W), got shape={img.shape}")

    rgb = to_uint8(rgb)
    rgb = np.ascontiguousarray(rgb)

    # Heuristic BigTIFF if output could exceed classic TIFF limits
    bigtiff = args.bigtiff or (rgb.nbytes > 2**31)

    tiff.imwrite(
        args.output_tif,
        rgb,
        photometric="rgb",
        compression=compression,
        bigtiff=bigtiff,
    )

    print(f"Wrote {args.output_tif}  shape={rgb.shape} dtype={rgb.dtype} bigtiff={bigtiff}")


if __name__ == "__main__":
    main()
