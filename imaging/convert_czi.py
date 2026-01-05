#!/usr/bin/env python
"""
Convert CZI files to OME-TIFF format
Usage: python convert_czi.py input.czi output.ome.tiff
"""

## conda create -n czi2tiff python=3.10 -y
## conda activate czi2tiff
## pip install aicsimageio[czi] ome-types

import sys
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter

def convert_czi_to_ometiff(input_path, output_path):
    """Convert a CZI file to OME-TIFF format"""
    print(f"Reading CZI file: {input_path}")
    img = AICSImage(input_path)
    
    print(f"Image dimensions: {img.dims}")
    print(f"Image shape: {img.shape}")
    print(f"Writing OME-TIFF: {output_path}")
    
    OmeTiffWriter.save(
        img.data,
        output_path,
        dim_order=img.dims.order,
        channel_names=img.channel_names,
        physical_pixel_sizes=img.physical_pixel_sizes
    )
    
    print("Conversion complete!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_czi.py input.czi output.ome.tiff")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_czi_to_ometiff(input_file, output_file)
