import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

# Input NetCDF file change accordingly
INPUT_FILE = r"C:\Users\lsh51\Downloads\Output\Output\GEOSChem_output\withoutfire(v1)\05x05_CEDS_1965_off_off_pm25_Surface_Re_yearavg.nc4"

# Output directory (#Desktop) change accordingly
DESKTOP = Path.home() / "Desktop"
OUTPUT_DIR = DESKTOP / "PM25_Data_2015"

OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("PM2.5 NetCDF Extraction Tool")
print("=" * 80)

print(f"\nInput file: {INPUT_FILE}")
print(f"Output directory: {OUTPUT_DIR}")

# =============================================================================
# Open NetCDF file
# =============================================================================

print("\n" + "─" * 80)
print("Opening NetCDF file...")
print("─" * 80)

try:
    dataset = nc.Dataset(INPUT_FILE, "r")
    print("✓ File successfully opened")
except FileNotFoundError:
    print("✗ Error: File not found")
    exit(1)
except Exception as e:
    print(f"✗ Error: {e}")
    exit(1)

# =============================================================================
# Save metadata
# =============================================================================

print("\n" + "─" * 80)
print("Saving metadata...")
print("─" * 80)

metadata_file = OUTPUT_DIR / "metadata.txt"

with open(metadata_file, "w", encoding="utf-8") as f:

    f.write("=" * 80 + "\n")
    f.write("PM2.5 NetCDF Metadata\n")
    f.write("=" * 80 + "\n\n")

    f.write(f"File name: {os.path.basename(INPUT_FILE)}\n\n")

    f.write("Global Attributes\n")
    f.write("-" * 80 + "\n")

    for attr in dataset.ncattrs():
        f.write(f"{attr}: {getattr(dataset, attr)}\n")

    f.write("\nDimensions\n")
    f.write("-" * 80 + "\n")

    for dim_name, dim in dataset.dimensions.items():
        f.write(f"{dim_name}: {len(dim)}\n")

    f.write("\nVariables\n")
    f.write("-" * 80 + "\n")

    for var_name, var in dataset.variables.items():

        f.write(f"\n{var_name}\n")
        f.write(f"  dimensions: {var.dimensions}\n")
        f.write(f"  shape: {var.shape}\n")
        f.write(f"  dtype: {var.dtype}\n")

print(f"✓ Metadata saved: {metadata_file}")

# =============================================================================
# Extract PM2.5
# =============================================================================

print("\n" + "─" * 80)
print("Extracting PM2.5 data...")
print("─" * 80)

pm25_var_name = None

for var in dataset.variables.keys():
    if "PM25" in var.upper() or "PM2.5" in var.upper():
        pm25_var_name = var
        break

if pm25_var_name is None:
    print("✗ Error: PM2.5 variable not found")
    print("Available variables:", list(dataset.variables.keys()))
    dataset.close()
    exit(1)

print(f"✓ PM2.5 variable found: {pm25_var_name}")

pm25 = dataset.variables[pm25_var_name][:]

print(f"Original PM2.5 shape: {pm25.shape}")

# Remove singleton dimensions (e.g. (1,1,360,720) -> (360,720))
pm25_2d = np.squeeze(pm25)

print(f"Final data shape: {pm25_2d.shape}")

# =============================================================================
# Read coordinates
# =============================================================================

lat = None
lon = None

for name in ["lat", "latitude"]:
    if name in dataset.variables:
        lat = dataset.variables[name][:]
        print(f"Latitude variable: {name} ({lat.shape})")
        break

for name in ["lon", "longitude"]:
    if name in dataset.variables:
        lon = dataset.variables[name][:]
        print(f"Longitude variable: {name} ({lon.shape})")
        break

if lat is None or lon is None:
    print("✗ Error: Latitude/Longitude not found")
    dataset.close()
    exit(1)

# =============================================================================
# Save matrix CSV
# =============================================================================

print("\n" + "─" * 80)
print("Saving matrix CSV...")
print("─" * 80)

matrix_csv = OUTPUT_DIR / "pm25_matrix.csv"

df_matrix = pd.DataFrame(
    pm25_2d,
    index=lat,
    columns=lon
)

df_matrix.index.name = "latitude"
df_matrix.columns.name = "longitude"

df_matrix.to_csv(matrix_csv)

print(f"✓ Matrix CSV saved: {matrix_csv}")

# =============================================================================
# Save long-format CSV
# =============================================================================

print("\n" + "─" * 80)
print("Saving long format CSV...")
print("─" * 80)

lon_grid, lat_grid = np.meshgrid(lon, lat)

df_long = pd.DataFrame({
    "latitude": lat_grid.flatten(),
    "longitude": lon_grid.flatten(),
    "pm25": pm25_2d.flatten()
})

long_csv = OUTPUT_DIR / "pm25_long.csv"

df_long.to_csv(long_csv, index=False)

print(f"✓ Long CSV saved: {long_csv}")

# =============================================================================
# Statistics
# =============================================================================

print("\n" + "─" * 80)
print("Computing statistics...")
print("─" * 80)

valid = pm25_2d.flatten()
valid = valid[~np.isnan(valid)]

stats_file = OUTPUT_DIR / "statistics.txt"

with open(stats_file, "w") as f:

    f.write("PM2.5 Statistics\n\n")

    f.write(f"Min: {valid.min():.4f}\n")
    f.write(f"Max: {valid.max():.4f}\n")
    f.write(f"Mean: {valid.mean():.4f}\n")
    f.write(f"Median: {np.median(valid):.4f}\n")
    f.write(f"Std: {valid.std():.4f}\n")

print(f"✓ Statistics saved: {stats_file}")

# =============================================================================
# Save NumPy file
# =============================================================================

npz_file = OUTPUT_DIR / "pm25_data.npz"

np.savez_compressed(
    npz_file,
    pm25=pm25_2d,
    latitude=lat,
    longitude=lon
)

print(f"✓ NumPy file saved: {npz_file}")

# =============================================================================
# Close file
# =============================================================================

dataset.close()

print("\n" + "=" * 80)
print("✓ Processing complete")
print("=" * 80)

print("\nOutput files:")

for f in OUTPUT_DIR.iterdir():
    size = f.stat().st_size / 1024
    print(f" • {f.name} ({size:.1f} KB)")