import sys
import math

# === PARAMETERS ===
R_CRIT = 10.0                 # critical radius for interaction
ANGLE_SHIFT = 100.0          # angular shift per helix period (degrees)
Z_SHIFT = 1.49              # vertical shift per helix period

ANGLE_MIN = -ANGLE_SHIFT/2 -180           # lower bound for normalized angle (deg)
ANGLE_MAX = -ANGLE_SHIFT/2 +180            # upper bound (exclusive)

# === FUNCTIONS ===

def normalize_angle_z(angle, delta_z):
    """
    Normalize (angle, delta_z) into the fundamental domain:
    - angle ∈ [ANGLE_MIN, ANGLE_MAX)
    - delta_z ∈ [0, Z_SHIFT)
    using periodicity: (angle + ANGLE_SHIFT, delta_z - Z_SHIFT)
    """
    # Apply spiral translation to bring z into [0, Z_SHIFT)
    while delta_z < 0:
        angle -= ANGLE_SHIFT
        delta_z += Z_SHIFT
    while delta_z >= Z_SHIFT:
        angle += ANGLE_SHIFT
        delta_z -= Z_SHIFT

    # Bring angle into [ANGLE_MIN, ANGLE_MAX)
    while angle >= ANGLE_MAX: angle -= 360
    while angle <  ANGLE_MIN: angle += 360

    return angle, delta_z

def apply_point_symmetry_if_needed(angle, delta_z, same_orientation):
    """
    If orientations are equal, apply point reflection symmetry:
    (angle, delta_z) → (–angle – ANGLE_SHIFT, Z_SHIFT – delta_z)
    if the mirrored angle lies within the more restricted angular domain.
    """
    if not same_orientation:
        return angle, delta_z  # symmetry not applicable

    mirrored_angle = -angle - ANGLE_SHIFT
    mirrored_z = Z_SHIFT - delta_z

    # Wrap mirrored angle into domain
    while mirrored_angle >= ANGLE_MAX: mirrored_angle -= 360
    while mirrored_angle  < ANGLE_MIN: mirrored_angle += 360

    # Accept mirrored values if within refined domain
    if ANGLE_MAX-180 <= mirrored_angle < ANGLE_MAX :
        return mirrored_angle, mirrored_z
    else:
        return angle, delta_z

# === MAIN PROGRAM ===

def main():
    if len(sys.argv) < 2:
        print("Error: Please provide a filename as a command-line argument.")
        sys.exit(1)

    filename = sys.argv[1]

    # Data containers
    x_list, y_list = [], []
    phi_list, z_list = [], []
    handedness_list, orientation_list = [], []

    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip()]
            if not lines:
                print("Error: File is empty.")
                return

            # Read unit cell dimensions
            unit_cell = list(map(float, lines[0].split()))
            if len(unit_cell) != 2:
                print("Error: First line must contain exactly 2 numbers for the unit cell.")
                return
            cell_x, cell_y = unit_cell
            print(f"Unit cell dimensions: x = {cell_x}, y = {cell_y}")

            # Read all entries
            for i, line in enumerate(lines[1:], start=2):
                parts = line.split()
                if len(parts) != 7:
                    print(f"Warning: Line {i} skipped — expected 7 values, got {len(parts)}")
                    continue

                x, y, phi, z, hand, orient, _ = map(float, parts)
                x_list.append(x)
                y_list.append(y)
                phi_list.append(phi)
                z_list.append(z)
                handedness_list.append(int(hand))
                orientation_list.append(int(orient))

        num_points = len(x_list)
        print(f"Successfully read {num_points} entries.\n")
        print("Pair interactions (i, j): Δz | angle | orientation code")

        # Pairwise interaction loop
        for i in range(num_points):
            for j in range(i + 1, num_points):
                dx = x_list[j] - x_list[i]
                dy = y_list[j] - y_list[i]
                distance = math.sqrt(dx**2 + dy**2)

                if distance < R_CRIT:
                    alpha_rad = math.atan2(dy, dx)
                    alpha_deg = math.degrees(alpha_rad)

                    phi_i_corr =  handedness_list[i]*orientation_list[i]*(phi_list[i] - alpha_deg)
                    phi_j_corr =  handedness_list[i]*orientation_list[i]*(phi_list[j] - alpha_deg)	
                    delta_z    = orientation_list[i]*(z_list[j] - z_list[i])    
                    
                    same_handed = handedness_list[i] == handedness_list[j]
                    same_orient = orientation_list[i] == orientation_list[j]

                    if same_handed:
                        angle_value = phi_i_corr - phi_j_corr
                        angle_type = "diff"
                    else:
                        angle_value = phi_j_corr + phi_i_corr
                        angle_type = "sum"
                        
                    

                    # Determine abbreviation
                    if same_handed and same_orient:
                        abbrev = "SS"
                    elif same_handed and not same_orient:
                        abbrev = "SD"
                    elif not same_handed and same_orient:
                        abbrev = "DS"
                    else:
                        abbrev = "DD"

                    # Normalize and apply point symmetry if needed
                    angle_norm, z_norm = normalize_angle_z(angle_value, delta_z)
                    angle_norm, z_norm = apply_point_symmetry_if_needed(angle_norm, z_norm, same_orient)

                    if angle_norm>ANGLE_MAX: print("now",file=sys.stderr)

                    # Output
                    print(f"{i:2d}, {j:2d} | dz = {z_norm:.3f} | "
                          f"{angle_type} = {angle_norm:.2f} | {abbrev} | {phi_i_corr%20:.3f}")

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()
