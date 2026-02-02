import tifffile as tiff
import sys

def extract_first_channel(input_file, output_file):
    try:
        # Open the TIFF file
        with tiff.TiffFile(input_file) as tif:
            # Access the first page corresponding to the first channel
            first_page = tif.pages[0]
            
            # Read only the first page's data into memory
            first_channel = first_page.asarray()

        # Save the first channel as a new TIFF
        tiff.imwrite(output_file, first_channel, photometric='minisblack')
        print(f"First channel saved to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_first_channel.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        extract_first_channel(input_file, output_file)
