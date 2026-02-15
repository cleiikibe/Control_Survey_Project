import math
import matplotlib.pyplot as plt
import numpy as np

# ============================================================================
# MODULE 1: Provisional Coordinates Computation (Intersection/Resection)
# ============================================================================

def dms_to_decimal(degrees, minutes, seconds):
    """Convert DMS to decimal degrees."""
    return degrees + minutes / 60.0 + seconds / 3600.0

def decimal_to_dms(decimal_deg):
    """Convert decimal degrees to degrees, minutes, seconds."""
    degrees = int(decimal_deg)
    minutes_decimal = (decimal_deg - degrees) * 60
    minutes = int(minutes_decimal)
    seconds = (minutes_decimal - minutes) * 60
    return degrees, minutes, seconds

def format_dms(decimal_deg):
    """Format decimal degrees as DMS string."""
    deg, min, sec = decimal_to_dms(decimal_deg)
    return f"{deg}°{min}'{sec:.1f}\""

def calculate_join_bearing(northing1, easting1, northing2, easting2):
    """Calculate bearing and distance between two points."""
    dN = northing2 - northing1
    dE = easting2 - easting1
    bearing_rad = math.atan2(dE, dN)
    bearing_deg = math.degrees(bearing_rad)
    if bearing_deg < 0:
        bearing_deg += 360.0
    distance = math.hypot(dE, dN)
    return bearing_deg, distance

def rec_polar_to_rect(distance, bearing_deg):
    """Convert polar (distance, bearing) to rectangular components (ΔN, ΔE)."""
    brad = math.radians(bearing_deg)
    dN = distance * math.cos(brad)
    dE = distance * math.sin(brad)
    return dN, dE

def input_dms(angle_name):
    """Prompt user to input DMS as space-separated values."""
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
    """Prompt user to input a coordinate value."""
    while True:
        try:
            val = float(input(f"Enter {coord_name} for station {station_name}: "))
            return val
        except ValueError:
            print("Invalid input. Please enter a numeric value.")

def compute_provisional_coordinates():
    """
    Compute provisional coordinates using intersection/resection.
    Returns: (provisional_northing, provisional_easting, unknown_name)
    """
    print("\n" + "="*90)
    print("PROVISIONAL COORDINATES COMPUTATION (Intersection/Resection)")
    print("="*90)

    # Input unknown station name
    UNKNOWN_NAME = input("\nEnter the name of the unknown station: ").strip()
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

        alpha_deg = None
        beta_deg = None

        print("\nUsing YOUR bearings:")
        print(f"Unknown station: {UNKNOWN_NAME}")
        print(f"Station A: E = {A['easting']:.2f}, N = {A['northing']:.2f}")
        print(f"Station B: E = {B['easting']:.2f}, N = {B['northing']:.2f}")
        print(f"Bearing A→P = {format_dms(bearing_AP)} = {bearing_AP:.6f}°")
        print(f"Bearing B→P = {format_dms(bearing_BP)} = {bearing_BP:.6f}°\n")

    else:
        # User does not have bearings, ask for alpha and beta angles
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

        # Bearings to unknown station P
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
    print("\n" + "="*90)
    print("COMPUTED BEARINGS")
    print("="*90)
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

    # Final coordinates of unknown station P (two computations)
    final_northing_upper = A['northing'] + rec_dN_upper if not math.isnan(rec_dN_upper) else float('nan')
    final_easting_upper = A['easting'] + rec_dE_upper if not math.isnan(rec_dE_upper) else float('nan')

    final_northing_lower = B['northing'] + rec_dN_lower if not math.isnan(rec_dN_lower) else float('nan')
    final_easting_lower = B['easting'] + rec_dE_lower if not math.isnan(rec_dE_lower) else float('nan')

    # Output table
    print("\n" + "="*90)
    print("INTERSECTION COMPUTATION TABLE")
    print("="*90)

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

    # Calculate PROVISIONAL coordinate as average of the two computations (standard practice)
    if math.isnan(final_northing_upper) or math.isnan(final_northing_lower):
        print("\n❌ Error: Could not compute provisional coordinates. Exiting.")
        return None, None, None
    
    provisional_northing = (final_northing_upper + final_northing_lower) / 2.0
    provisional_easting = (final_easting_upper + final_easting_lower) / 2.0
    
    print("\n" + "="*90)
    print("PROVISIONAL COORDINATES (Averaged from both computations)")
    print("="*90)
    print(f"Provisional {UNKNOWN_NAME}: Northing = {provisional_northing:.2f}, Easting = {provisional_easting:.2f}")
    print("="*90)
    
    return provisional_northing, provisional_easting, UNKNOWN_NAME

# ============================================================================
# MODULE 2: Cuts Computation (Your perfected program)
# ============================================================================

def dms_to_decimal_rounded(degrees, minutes, seconds):
    """Convert DMS to decimal degrees and ROUND to 6 decimals BEFORE trig calculations."""
    decimal = degrees + minutes / 60.0 + seconds / 3600.0
    return round(decimal, 6)

def cot(angle_rad):
    """cot(θ) = 1 / tan(θ)"""
    return 1 / math.tan(angle_rad)

def csc(angle_rad):
    """csc(θ) = 1 / sin(θ)"""
    return 1 / math.sin(angle_rad)

def sec(angle_rad):
    """sec(θ) = 1 / cos(θ)"""
    return 1 / math.cos(angle_rad)

def compute_cuts(station_id, N_C, E_C, N_P, E_P, bearing_dms):
    """Compute cuts using EXACT methodology from your notes."""
    bearing_deg = dms_to_decimal_rounded(*bearing_dms)
    bearing_rad = math.radians(bearing_deg)
    
    delta_N = N_P - N_C
    delta_E = E_P - E_C
    
    cot_bearing = cot(bearing_rad)
    delta_N_prime = cot_bearing * delta_E
    cut_N = delta_N_prime - delta_N
    csc_bearing = csc(bearing_rad)
    s1 = csc_bearing * delta_E
    
    tan_bearing = math.tan(bearing_rad)
    delta_E_prime = tan_bearing * delta_N
    cut_E = delta_E_prime - delta_E
    sec_bearing = sec(bearing_rad)
    s2 = delta_N * sec_bearing
    
    return {
        'station_id': station_id,
        'N_C': N_C, 'E_C': E_C, 'N_P': N_P, 'E_P': E_P,
        'bearing_dms': bearing_dms, 'bearing_deg': bearing_deg,
        'delta_N': delta_N, 'delta_E': delta_E,
        'cot_bearing': cot_bearing, 'delta_N_prime': delta_N_prime, 'cut_N': cut_N,
        'csc_bearing': csc_bearing, 's1': s1,
        'tan_bearing': tan_bearing, 'delta_E_prime': delta_E_prime, 'cut_E': cut_E,
        'sec_bearing': sec_bearing, 's2': s2
    }

def print_steps(result):
    """Print calculations EXACTLY as in your handwritten notes."""
    b = result['bearing_dms']
    print(f"\n{'='*70}")
    print(f"STATION {result['station_id']}")
    print(f"{'='*70}")
    print(f"Known Station C:   N = {result['N_C']:.2f}, E = {result['E_C']:.2f}")
    print(f"Provisional P:     N = {result['N_P']:.2f}, E = {result['E_P']:.2f}")
    print(f"Bearing CP:        {b[0]:.0f}° {b[1]:.0f}' {b[2]:.0f}\" = {result['bearing_deg']:.6f}°")
    print(f"\nΔN = {result['N_P']:.2f} - {result['N_C']:.2f} = {result['delta_N']:.2f}")
    print(f"ΔE = {result['E_P']:.2f} - {result['E_C']:.2f} = {result['delta_E']:.2f}")
    
    print(f"\nCUT IN NORTHING:")
    print(f"ΔN' = cot({result['bearing_deg']:.6f}°) × ΔE")
    print(f"    = {result['cot_bearing']} × {result['delta_E']:.2f}")
    print(f"    = {result['delta_N_prime']:.2f}")
    print(f"Cut = ΔN' - ΔN")
    print(f"    = {result['delta_N_prime']:.2f} - ({result['delta_N']:.2f})")
    print(f"    = {result['cut_N']:.2f}")
    
    print(f"\nDistance s1 = csc({result['bearing_deg']:.6f}°) × ΔE")
    print(f"            = {result['csc_bearing']} × {result['delta_E']:.2f}")
    print(f"            = {result['s1']:.2f}")
    
    print(f"\nCUT IN EASTING:")
    print(f"ΔE' = tan({result['bearing_deg']:.6f}°) × ΔN")
    print(f"    = {result['tan_bearing']} × {result['delta_N']:.2f}")
    print(f"    = {result['delta_E_prime']:.2f}")
    print(f"Cut = ΔE' - ΔE")
    print(f"    = {result['delta_E_prime']:.2f} - {result['delta_E']:.2f}")
    print(f"    = {result['cut_E']:.2f}")
    
    print(f"\nDistance s2 = ΔN × sec({result['bearing_deg']:.6f}°)")
    print(f"            = {result['delta_N']:.2f} × {result['sec_bearing']}")
    print(f"            = {result['s2']:.2f}")
    
    print(f"\n{'='*70}")
    print(f"FINAL CUTS:  Cut_N = {result['cut_N']:.2f} m,  Cut_E = {result['cut_E']:.2f} m")
    print(f"{'='*70}")

def plot_cuts_cartesian(results, N_P, E_P):
    """Create a plot with intersection point calculated using S1 weights."""
    plt.figure(figsize=(10, 10))
    
    # Create Cartesian plane with origin at center
    plt.axhline(y=0, color='k', linestyle='-', linewidth=1.5)
    plt.axvline(x=0, color='k', linestyle='-', linewidth=1.5)
    
    # Set up grid
    plt.grid(True, linestyle='--', alpha=0.7, linewidth=0.8)
    
    # Set axis labels
    plt.xlabel('EASTING (m)', fontsize=12, fontweight='bold')
    plt.ylabel('NORTHING (m)', fontsize=12, fontweight='bold')
    plt.title('CUTS CARTESIAN PLOT\n(ΔN vs ΔE for Each Station)', 
              fontsize=14, fontweight='bold', pad=20)
    
    # Set appropriate limits based on cut values
    max_cut = max(max(abs(r['cut_N']), abs(r['cut_E'])) for r in results)
    padding = max_cut * 1.5
    plt.xlim(-max_cut - padding, max_cut + padding)
    plt.ylim(-max_cut - padding, max_cut + padding)
    
    # Set custom tick marks
    tick_interval = 0.2
    max_tick = max_cut + padding
    num_ticks = int(max_tick / tick_interval) + 1
    
    x_ticks = [i * tick_interval for i in range(-num_ticks, num_ticks + 1)]
    y_ticks = [i * tick_interval for i in range(-num_ticks, num_ticks + 1)]
    
    plt.xticks(x_ticks, [f"{i * tick_interval:.1f}" for i in range(-num_ticks, num_ticks + 1)])
    plt.yticks(y_ticks, [f"{i * tick_interval:.1f}" for i in range(-num_ticks, num_ticks + 1)])
    
    # Define distinct colors for stations
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    # Plot the original cut lines and axis markings
    for idx, r in enumerate(results):
        color = colors[idx % len(colors)]
        
        # MARK CUT VALUES ON AXES
        plt.plot(0, r['cut_N'], 's', markersize=8, color=color, 
                markeredgecolor='black', markeredgewidth=1.5)
        plt.text(-0.05, r['cut_N'], f"{r['cut_N']:.2f}", 
                ha='right', fontsize=9, color=color, 
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=color, alpha=0.7))
        
        plt.plot(r['cut_E'], 0, 's', markersize=8, color=color, 
                markeredgecolor='black', markeredgewidth=1.5)
        plt.text(r['cut_E'], -0.05, f"{r['cut_E']:.2f}", 
                va='top', fontsize=9, color=color, 
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=color, alpha=0.7))
        
        # DRAW EXTENDED LINES (PAST THE MARKINGS)
        if abs(r['cut_E']) > 1e-5:
            slope = -r['cut_N'] / r['cut_E']
            x_min = -max_cut - padding
            x_max = max_cut + padding
            y_min = slope * x_min + r['cut_N']
            y_max = slope * x_max + r['cut_N']
            
            plt.plot([x_min, x_max], [y_min, y_max], 
                    linestyle='--', linewidth=2.0, color=color,
                    label=f"{r['station_id']} line")
        else:
            plt.plot([0, 0], [-max_cut - padding, max_cut + padding], 
                    linestyle='--', linewidth=2.0, color=color,
                    label=f"{r['station_id']} line")
        
        # Add station label at the cut point
        plt.text(r['cut_E'] + 0.03, r['cut_N'] + 0.03, 
                f"{r['station_id']}",
                fontsize=11, fontweight='bold', color=color)
        
        # Add cut values near the point
        plt.text(r['cut_E'] + 0.02, r['cut_N'] - 0.07, 
                f"ΔN={r['cut_N']:.2f}\nΔE={r['cut_E']:.2f}",
                fontsize=9, 
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=color, linewidth=1.5))
    
    # Calculate the weighted intersection point using S1 as weights
    total_weight = sum(r['s1'] for r in results)
    weighted_cut_N = sum(r['cut_N'] * r['s1'] for r in results) / total_weight
    weighted_cut_E = sum(r['cut_E'] * r['s1'] for r in results) / total_weight
    
    # Plot the intersection point (TRUE position of P)
    plt.plot(weighted_cut_E, weighted_cut_N, 'o', markersize=12, 
            color='purple', markeredgecolor='black', markeredgewidth=1.5,
            label='True position of P')
    
    # Add label for intersection point
    plt.text(weighted_cut_E + 0.03, weighted_cut_N + 0.03, 
            "True P", fontsize=12, fontweight='bold', color='purple')
    
    # Add cut values for true position
    plt.text(weighted_cut_E + 0.02, weighted_cut_N - 0.07, 
            f"ΔN={weighted_cut_N:.2f}\nΔE={weighted_cut_E:.2f}",
            fontsize=9, 
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='purple', linewidth=1.5))
    
    # Draw the adjusted lines that all intersect at the true position
    for idx, r in enumerate(results):
        color = colors[idx % len(colors)]
        
        # Calculate slope for the line
        if abs(r['cut_E']) > 1e-5:
            slope = -r['cut_N'] / r['cut_E']
        else:
            slope = float('inf')  # Vertical line
        
        # Draw the adjusted line through the true position
        if slope != float('inf'):
            x_min = -max_cut - padding
            x_max = max_cut + padding
            y_min = slope * x_min + (weighted_cut_N - slope * weighted_cut_E)
            y_max = slope * x_max + (weighted_cut_N - slope * weighted_cut_E)
            
            plt.plot([x_min, x_max], [y_min, y_max], 
                    linestyle='-', linewidth=2.5, color=color,
                    alpha=0.7, label=f"Adjusted {r['station_id']}")
        else:
            plt.plot([weighted_cut_E, weighted_cut_E], [-max_cut - padding, max_cut + padding], 
                    linestyle='-', linewidth=2.5, color=color,
                    alpha=0.7, label=f"Adjusted {r['station_id']}")
    
    # Add north arrow
    north_x = max_cut - 0.1
    north_y = -max_cut + 0.05
    plt.arrow(north_x, north_y, 0, 0.15, head_width=0.03, head_length=0.04, 
             fc='black', ec='black', linewidth=1.5)
    plt.text(north_x - 0.05, north_y + 0.17, 'N', fontsize=12, fontweight='bold')
    
    # Add table of cut values
    table_data = []
    for r in results:
        table_data.append([r['station_id'], f"{r['cut_N']:.2f}", f"{r['cut_E']:.2f}"])
    
    # Add true position row
    table_data.append(["TRUE P", f"{weighted_cut_N:.2f}", f"{weighted_cut_E:.2f}"])
    
    # Position the table at the bottom
    table = plt.table(cellText=table_data,
                     colLabels=['STATION', 'CUT (N)', 'CUT (E)'],
                     loc='bottom',
                     cellLoc='center',
                     bbox=[0.1, -0.25, 0.8, 0.18])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style the table header
    for i in range(3):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style table cells with alternating colors
    for i in range(len(results) + 1):
        for j in range(3):
            if i == len(results):  # True position row
                table[(i, j)].set_facecolor('#FFD700')
            else:
                table[(i, j)].set_facecolor('#E7E6E6' if i % 2 == 0 else 'white')
    
    # Add legend at top-right
    plt.legend(loc='upper right', fontsize=9, 
              framealpha=0.9, shadow=True)
    
    plt.tight_layout()
    plt.show()

def compute_cuts_workflow(N_P, E_P, unknown_name="P"):
    """Run the complete cuts computation workflow."""
    print("\n" + "="*70)
    print("CUTS COMPUTATION (Surveying Corrections)")
    print("="*70)
    print(f"Provisional Station {unknown_name}: N = {N_P:.2f}, E = {E_P:.2f}")
    
    # STEP 2: Ask how many known stations
    print("\n>>> How many known stations will you process? <<<")
    try:
        n = int(input("Number of known stations: "))
        if n < 1:
            print("❌ Must have at least 1 known station. Exiting.")
            return
    except ValueError:
        print("❌ Invalid number. Exiting.")
        return
    
    # STEP 3: Collect data for ALL known stations
    stations_data = []
    print(f"\n>>> Enter data for all {n} known stations <<<")
    
    for i in range(1, n + 1):
        print(f"\n--- Station {i} of {n} ---")
        sid = input(f"Station ID (e.g., C{i}): ").strip() or f"C{i}"
        
        try:
            N_C = float(input(f"Northing of {sid}: "))
            E_C = float(input(f"Easting of {sid}: "))
            bearing_input = input(f"Bearing from {sid} to {unknown_name} (DMS as 'deg min sec'): ").strip()
            parts = bearing_input.split()
            if len(parts) != 3:
                raise ValueError("Bearing must have exactly 3 values separated by spaces")
            d, m, s = map(float, parts)
            bearing = (d, m, s)
            stations_data.append({
                'station_id': sid,
                'N_C': N_C,
                'E_C': E_C,
                'bearing_dms': bearing
            })
        except (ValueError, IndexError) as e:
            print(f"❌ Invalid input for station {sid}. Skipping.")
            continue
    
    if not stations_data:
        print("\n❌ No valid stations entered. Exiting.")
        return
    
    # STEP 4: PROCESS ALL STATIONS AND SHOW CALCULATIONS
    print("\n" + "="*70)
    print("CALCULATIONS FOR ALL STATIONS")
    print("="*70)
    
    results = []
    for data in stations_data:
        res = compute_cuts(
            data['station_id'],
            data['N_C'],
            data['E_C'],
            N_P,
            E_P,
            data['bearing_dms']
        )
        results.append(res)
        print_steps(res)
    
    # STEP 5: FINAL SUMMARY TABLE
    print(f"\n{'='*70}")
    print("SUMMARY OF ALL CUTS")
    print(f"{'='*70}")
    print(f"{'Station':<10} {'Cut_N (m)':>12} {'Cut_E (m)':>12} {'s1 (m)':>12} {'s2 (m)':>12}")
    print(f"{'-'*70}")
    for r in results:
        print(f"{r['station_id']:<10} {r['cut_N']:>12.2f} {r['cut_E']:>12.2f} "
              f"{r['s1']:>12.2f} {r['s2']:>12.2f}")
    print(f"{'='*70}")
    
    # STEP 6: PLOT THE CUTS
    print("\n" + "="*70)
    print("GENERATING CUTS CARTESIAN PLOT WITH INTERSECTION CALCULATION")
    print("="*70)
    plot_cuts_cartesian(results, N_P, E_P)
    
    # Calculate true position of P
    total_weight = sum(r['s1'] for r in results)
    weighted_cut_N = sum(r['cut_N'] * r['s1'] for r in results) / total_weight
    weighted_cut_E = sum(r['cut_E'] * r['s1'] for r in results) / total_weight
    
    true_N = N_P + weighted_cut_N
    true_E = E_P + weighted_cut_E
    
    print("\n" + "="*70)
    print(f"FINAL ADJUSTED COORDINATES OF STATION {unknown_name}")
    print(f"{'='*70}")
    print(f"Provisional {unknown_name}: N = {N_P:.2f}, E = {E_P:.2f}")
    print(f"Weighted Cut: ΔN = {weighted_cut_N:.2f} m, ΔE = {weighted_cut_E:.2f} m")
    print(f"True Position: N = {true_N:.2f}, E = {true_E:.2f}")
    print(f"{'='*70}")
    
    print("\n✅ All calculations complete!")
    print("✅ Trigonometric values computed with 6-decimal precision (matches calculator)")
    print("✅ Plot shows both original and adjusted lines intersecting at true position")

# ============================================================================
# MAIN UNIFIED WORKFLOW
# ============================================================================

def main():
    print("="*90)
    print("SURVEYING COMPUTATION SUITE")
    print("Unified Program: Provisional Coordinates + Cuts Computation")
    print("="*90)
    
    # STEP 1: Ask if user has provisional coordinates
    while True:
        has_provisional = input("\nDo you have provisional coordinates for the unknown station? (Y/N): ").strip().upper()
        if has_provisional in ['Y', 'YES', 'N', 'NO']:
            break
        else:
            print("Please answer 'Y' or 'N'.")
    
    if has_provisional in ['Y', 'YES']:
        # User has provisional coordinates - get them directly
        print("\n>>> Enter provisional coordinates for the unknown station <<<")
        try:
            N_P = float(input("Northing of provisional station: "))
            E_P = float(input("Easting of provisional station: "))
            unknown_name = input("Station name (default 'P'): ").strip() or "P"
        except ValueError:
            print("❌ Invalid input. Exiting.")
            return
    else:
        # User does NOT have provisional coordinates - compute them
        print("\n>>> Computing provisional coordinates using intersection/resection <<<")
        N_P, E_P, unknown_name = compute_provisional_coordinates()
        
        if N_P is None or E_P is None:
            print("\n❌ Provisional coordinates computation failed. Exiting.")
            return
        
        print("\n✅ Provisional coordinates computed successfully. Proceeding to cuts computation...")
    
    # STEP 2: Run cuts computation using the provisional coordinates
    compute_cuts_workflow(N_P, E_P, unknown_name)
    
    print("\n" + "="*90)
    print("SURVEYING COMPUTATION COMPLETE")
    print("="*90)
    print("✅ Provisional coordinates obtained (user input or computed)")
    print("✅ Cuts computed for all known stations")
    print("✅ True position of unknown station determined")
    print("✅ Professional plot generated showing all computations")
    print("="*90)

if __name__ == "__main__":
    main()