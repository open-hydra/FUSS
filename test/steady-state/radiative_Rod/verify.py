import os
import re
import sys

WALL_TEC     = os.path.join(os.path.dirname(__file__), "OUTPUT", "wall.tec")
TARGET_ZONE  = "B1F2"
REFERENCE_TW = 927.0 # K


def parse_wall_tec(filepath):
    
    with open(filepath, "r") as fh:
        lines = fh.readlines()

    variables = []
    zones = {}
    current_header = None
    raw = []

    for line in lines:
        stripped = line.strip()

        if re.match(r"(?i)variables", stripped):
            variables = re.findall(r'"([^"]+)"', stripped)
            continue

        if re.match(r"(?i)zone", stripped):
            if current_header is not None:
                zones[current_header["name"]] = {"header": current_header, "raw": raw}
            raw = []

            name_m = re.search(r'T\s*=\s*(\w+)',   stripped, re.IGNORECASE)
            i_m    = re.search(r'\bI\s*=\s*(\d+)', stripped, re.IGNORECASE)
            j_m    = re.search(r'\bJ\s*=\s*(\d+)', stripped, re.IGNORECASE)
            k_m    = re.search(r'\bK\s*=\s*(\d+)', stripped, re.IGNORECASE)
            cc_m   = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED',
                               stripped, re.IGNORECASE)

            current_header = {
                "name"    : name_m.group(1) if name_m else f"Zone{len(zones)+1}",
                "ni"      : int(i_m.group(1)) if i_m else 1,
                "nj"      : int(j_m.group(1)) if j_m else 1,
                "nk"      : int(k_m.group(1)) if k_m else 1,
                "cc_start": int(cc_m.group(1)) if cc_m else 0,
                "cc_end"  : int(cc_m.group(2)) if cc_m else 0,
            }
            continue

        try:
            raw.append(float(stripped))
        except ValueError:
            pass

    if current_header is not None:
        zones[current_header["name"]] = {"header": current_header, "raw": raw}

    return variables, zones


def extract_vars(variables, zone_info):
    
    hdr = zone_info["header"]
    raw = zone_info["raw"]
    ni, nj, nk       = hdr["ni"], hdr["nj"], hdr["nk"]
    cc_start, cc_end = hdr["cc_start"], hdr["cc_end"]

    n_nodal = ni * nj * nk
    n_cell  = max(ni-1,1) * max(nj-1,1) * max(nk-1,1)

    result = {}
    idx = 0
    for v_idx, vname in enumerate(variables, start=1):
        n = n_cell if (cc_start <= v_idx <= cc_end) else n_nodal
        result[vname] = raw[idx : idx + n]
        idx += n

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    if not os.path.isfile(WALL_TEC):
        print(f"ERROR: file not found – {WALL_TEC}")
        sys.exit(1)

    variables, zones = parse_wall_tec(WALL_TEC)

    if TARGET_ZONE not in zones:
        available = ", ".join(zones.keys())
        print(f"ERROR: zone '{TARGET_ZONE}' not found. Available zones: {available}")
        sys.exit(1)

    zone_vars = extract_vars(variables, zones[TARGET_ZONE])

    if "Tw" not in zone_vars:
        print(f"ERROR: variable 'Tw' not found in zone {TARGET_ZONE}.")
        print(f"       Available variables: {', '.join(variables)}")
        sys.exit(1)

    tw_values  = zone_vars["Tw"]
    tw_fuss    = sum(tw_values) / len(tw_values)
    error_abs  = tw_fuss - REFERENCE_TW
    error_pct  = abs(error_abs) / REFERENCE_TW * 100.0

    sep = "-" * 48

    print()
    print(sep)
    print(" 1-D STEADY STATE THERMAL ANALYSIS ")
    print(" NAFEMS Benchmark T9 (iii) ")
    print(sep)
    print(f"  {'Quantity':<20} {'Value (K)':>10}  {'Value (°C)':>10}")
    print(sep)
    print(f"  {'NAFEMS Reference':<20} "
          f"{REFERENCE_TW:>10.3f}  {REFERENCE_TW - 273.15:>10.3f}")
    print(f"  {'FUSS Solution':<20} "
          f"{tw_fuss:>10.3f}  {tw_fuss - 273.15:>10.3f}")
    print(sep)
    print(f"  {'Absolute error':<20} {error_abs:>10.3f}  {'K':>10}")
    print(f"  {'Relative error':<20} {error_pct:>10.4f}  {'%':>10}")
    print(sep)
    print()

    if error_pct < 1.0:
        print("  Result: PASS  (relative error < 1 %)")
    else:
        print("  Result: FAIL  (relative error >= 1 %)")

    print()


if __name__ == "__main__":
    main()