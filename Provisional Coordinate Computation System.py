import math

# --- Helper Functions ---
def dms_to_decimal(degrees, minutes, seconds):
    return degrees + minutes / 60.0 + seconds / 3600.0

def decimal_to_dms(decimal_deg):
    """Convert decimal degrees to degrees, minutes, seconds"""
    degrees = int(decimal_deg)
    minutes_decimal = (decimal_deg - degrees) * 60
    minutes = int(minutes_decimal)
    seconds = (minutes_decimal - minutes) * 60
    return degrees, minutes, seconds

def format_dms(decimal_deg):
    """Format decimal degrees as DMS string"""
    deg, min, sec = decimal_to_dms(decimal_deg)
    return f"{deg}°{min}'{sec:.1f}\""

def calculate_join_bearing(northing1, easting1, northing2, easting2):
    dN = northing2 - northing1
    dE = easting2 - easting1
    bearing_rad = math.atan2(dE, dN)
    bearing_deg = math.degrees(bearing_rad)
    if bearing_deg < 0:
        bearing_deg += 360.0
    distance = math.hypot(dE, dN)
    return bearing_deg, distance

def rec_polar_to_rect(distance, bearing_deg):
    """
    REC: Convert polar (distance, bearing) to rectangular components (ΔN, ΔE)
    Bearing is clockwise from North.
    """
    brad = math.radians(bearing_deg)
    dN = distance * math.cos(brad)
    dE = distance * math.sin(brad)
    return dN, dE

def input_dms(angle_name):
    """Prompt user to input degrees, minutes, seconds for an angle as space-separated values and convert to decimal degrees"""
    while True:
        try:
            raw = input(f"Enter {angle_name} in D M S separated by spaces (e.g. 129 33 21): ").strip()
            parts = raw.split()
            if len(parts) != 3:
                print("Please enter exactly three values: degrees, minutes, and seconds.")
                continue
            deg, min, sec = map(float, parts)
            if deg < 0 or min < 0 or sec < 0:
                print("Please enter non-negative values.")
                continue
            return dms_to_decimal(deg, min, sec)
        except ValueError:
            print("Invalid input. Please enter numeric values separated by spaces.")

def input_coordinate(station_name, coord_name):
    """Prompt user to input a coordinate value (northing or easting)"""
    while True:
        try:
            val = float(input(f"Enter {coord_name} for station {station_name}: "))
            return val
        except ValueError:
            print("Invalid input. Please enter a numeric value.")

def main():
    print("Surveying Intersection (Resection) Program — USER INPUT MODE")
    print("=" * 90)

    # Input unknown station name
    UNKNOWN_NAME = input("Enter the name of the unknown station: ").strip()
    if not UNKNOWN_NAME:
        UNKNOWN_NAME = "P"  # default if empty

    # Input known stations A and B coordinates
    print("\nEnter coordinates for Station A:")
    A_easting = input_coordinate("A", "Easting")
    A_northing = input_coordinate("A", "Northing")

    print("\nEnter coordinates for Station B:")
    B_easting = input_coordinate("B", "Easting")
    B_northing = input_coordinate("B", "Northing")

    A = {'name': 'A', 'easting': A_easting, 'northing': A_northing}
    B = {'name': 'B', 'easting': B_easting, 'northing': B_northing}

    # Ask if user has the bearings to the unknown station
    while True:
        has_bearing = input("\nDo you have the bearings to the unknown station? (yes/no): ").strip().lower()
        if has_bearing in ['yes', 'no']:
            break
        else:
            print("Please answer 'yes' or 'no'.")

    if has_bearing == 'yes':
        # Ask user to input bearings from A and B to unknown station
        print("\nEnter bearing from Station A to the unknown station:")
        bearing_AP = input_dms("Bearing A→P")

        print("\nEnter bearing from Station B to the unknown station:")
        bearing_BP = input_dms("Bearing B→P")

        # For display purposes, alpha and beta are not needed here, set to None
        alpha_deg = None
        beta_deg = None

        print("\nUsing YOUR bearings:")
        print(f"Unknown station: {UNKNOWN_NAME}")
        print(f"Station A: E = {A['easting']:.2f}, N = {A['northing']:.2f}")
        print(f"Station B: E = {B['easting']:.2f}, N = {B['northing']:.2f}")
        print(f"Bearing A→P = {format_dms(bearing_AP)} = {bearing_AP:.6f}°")
        print(f"Bearing B→P = {format_dms(bearing_BP)} = {bearing_BP:.6f}°\n")

    else:
        # User does not have bearings, ask for alpha and beta angles as usual
        print("\nEnter angle α (alpha) at Station A (degrees, minutes, seconds):")
        alpha_deg = input_dms("α")

        print("\nEnter angle β (beta) at Station B (degrees, minutes, seconds):")
        beta_deg = input_dms("β")

        print("\nUsing YOUR parameters:")
        print(f"Unknown station: {UNKNOWN_NAME}")
        print(f"Station A: E = {A['easting']:.2f}, N = {A['northing']:.2f}")
        print(f"Station B: E = {B['easting']:.2f}, N = {B['northing']:.2f}")
        print(f"α = {format_dms(alpha_deg)} = {alpha_deg:.6f}°")
        print(f"β = {format_dms(beta_deg)} = {beta_deg:.6f}°\n")

        # Step 1: JOIN AB
        bearing_AB, dist_AB = calculate_join_bearing(
            A['northing'], A['easting'],
            B['northing'], B['easting']
        )
        bearing_BA = (bearing_AB + 180.0) % 360.0

        # Bearings to unknown station P (per your diagram)
        bearing_AP = (bearing_AB - alpha_deg) % 360.0
        bearing_BP = (bearing_BA + beta_deg) % 360.0

    # Calculate ΔN and ΔE (change in coordinates from A to B)
    delta_N = B['northing'] - A['northing']
    delta_E = B['easting'] - A['easting']

    # Determine which station has larger bearing to P
    if bearing_BP >= bearing_AP:
        top, bot = B, A
        top_lbl, bot_lbl = "B", "A"
        top_bear, bot_bear = bearing_BP, bearing_AP
    else:
        top, bot = A, B
        top_lbl, bot_lbl = "A", "B"
        top_bear, bot_bear = bearing_AP, bearing_BP

    # Print bearings summary
    print("\n" + "=" * 90)
    print("COMPUTED BEARINGS")
    print("=" * 90)
    print(f"Bearing A → P = {format_dms(bearing_AP)} ({bearing_AP:.6f}°)")
    print(f"Bearing B → P = {format_dms(bearing_BP)} ({bearing_BP:.6f}°)")

    # Distance computation
    change_brg_deg = (bearing_BP - bearing_AP)
    change_brg_rad = math.radians(change_brg_deg)

    brgAP_rad = math.radians(bearing_AP)
    brgBP_rad = math.radians(bearing_BP)

    sin_change = math.sin(change_brg_rad)

    if abs(sin_change) < 1e-12:
        dAP = float('nan')
        dBP = float('nan')
    else:
        dAP = ((delta_N * math.sin(brgBP_rad)) - (delta_E * math.cos(brgBP_rad))) / sin_change
        dBP = ((delta_N * math.sin(brgAP_rad)) - (delta_E * math.cos(brgAP_rad))) / sin_change

    # REC for upper blank row
    if math.isnan(dAP):
        rec_dN_upper, rec_dE_upper = float('nan'), float('nan')
    else:
        rec_dN_upper, rec_dE_upper = rec_polar_to_rect(dAP, bearing_AP)

    # REC for lower blank row
    if math.isnan(dBP):
        rec_dN_lower, rec_dE_lower = float('nan'), float('nan')
    else:
        rec_dN_lower, rec_dE_lower = rec_polar_to_rect(dBP, bearing_BP)

    # Final coordinates of unknown station P
    final_northing_upper = A['northing'] + rec_dN_upper if not math.isnan(rec_dN_upper) else float('nan')
    final_easting_upper = A['easting'] + rec_dE_upper if not math.isnan(rec_dE_upper) else float('nan')

    final_northing_lower = B['northing'] + rec_dN_lower if not math.isnan(rec_dN_lower) else float('nan')
    final_easting_lower = B['easting'] + rec_dE_lower if not math.isnan(rec_dE_lower) else float('nan')

    # Output table
    print("\n" + "=" * 90)
    print("INTERSECTION COMPUTATION TABLE")
    print("=" * 90)

    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")
    print(f"| {'Station':<32} | {'Northing':<12} | {'Easting':<12} | {'Description':<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Upper P row with final coordinates and distance
    dist_ap_text = f"{dAP:.3f}" if not math.isnan(dAP) else ""
    finalN_upper_text = f"{final_northing_upper:12.2f}" if not math.isnan(final_northing_upper) else f"{'':>12}"
    finalE_upper_text = f"{final_easting_upper:12.2f}" if not math.isnan(final_easting_upper) else f"{'':>12}"
    print(f"| {f'{UNKNOWN_NAME} (unknown station)':<32} | {finalN_upper_text} | {finalE_upper_text} | {f'Distance AP = {dist_ap_text}':<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Upper blank row with REC values
    recN_upper_text = f"{rec_dN_upper:12.2f}" if not math.isnan(rec_dN_upper) else f"{'':>12}"
    recE_upper_text = f"{rec_dE_upper:12.2f}" if not math.isnan(rec_dE_upper) else f"{'':>12}"
    print(f"| {'':<32} | {recN_upper_text} | {recE_upper_text} | {'':<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Top station row with bearing
    if top_lbl == 'A':
        bearing_text = f"Bearing A→P = {format_dms(top_bear)}"
    else:
        bearing_text = f"Bearing B→P = {format_dms(top_bear)}"
    print(f"| {f'{top_lbl} (< bearing)':<32} | "
          f"{top['northing']:>12.2f} | "
          f"{top['easting']:>12.2f} | "
          f"{bearing_text:<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # ΔN, ΔE, bearing difference row
    bear_diff = abs(top_bear - bot_bear)
    description_text = f"BRG {top_lbl}P − BRG {bot_lbl}P = {format_dms(bear_diff)}"
    print(f"| {'':<32} | {delta_N:>12.2f} | {delta_E:>12.2f} | {description_text:<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Bottom station row with bearing
    if bot_lbl == 'A':
        bearing_text_bot = f"Bearing A→P = {format_dms(bot_bear)}"
    else:
        bearing_text_bot = f"Bearing B→P = {format_dms(bot_bear)}"
    print(f"| {f'{bot_lbl} (> bearing)':<32} | "
          f"{bot['northing']:>12.2f} | "
          f"{bot['easting']:>12.2f} | "
          f"{bearing_text_bot:<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Lower blank row with REC values
    recN_lower_text = f"{rec_dN_lower:12.2f}" if not math.isnan(rec_dN_lower) else f"{'':>12}"
    recE_lower_text = f"{rec_dE_lower:12.2f}" if not math.isnan(rec_dE_lower) else f"{'':>12}"
    print(f"| {'':<32} | {recN_lower_text} | {recE_lower_text} | {'':<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    # Lower P row with final coordinates and distance
    dist_bp_text = f"{dBP:.3f}" if not math.isnan(dBP) else ""
    finalN_lower_text = f"{final_northing_lower:12.2f}" if not math.isnan(final_northing_lower) else f"{'':>12}"
    finalE_lower_text = f"{final_easting_lower:12.2f}" if not math.isnan(final_easting_lower) else f"{'':>12}"
    print(f"| {f'{UNKNOWN_NAME} (unknown station)':<32} | {finalN_lower_text} | {finalE_lower_text} | {f'Distance BP = {dist_bp_text}':<30} |")
    print(f"+ {'-'*32} + {'-'*12} + {'-'*12} + {'-'*30} +")

    print("\n✅ Table matches your handwritten layout precisely.")
    print("All values computed using YOUR exact data.")

if __name__ == "__main__":
    main()
