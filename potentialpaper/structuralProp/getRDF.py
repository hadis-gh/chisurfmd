import sys
import numpy as np

def init_mode(rmin, rmax, dr):
    r_values = np.arange(rmin, rmax, dr)
    nbins = len(r_values)
    
    rdf = {
        "SS": np.zeros(nbins),
        "SD": np.zeros(nbins),
        "DS": np.zeros(nbins),
        "DD": np.zeros(nbins),
    }
    N = 0
    return N, r_values, rdf

def read_mode(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            print("Error: File does not contain enough data.")
            return 0, [], {}

        N = int(lines[0].strip())

        r_values = []
        rdf = {
            "SS": [],
            "SD": [],
            "DS": [],
            "DD": [],
        }

        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) != 5:
                continue
            r = float(parts[0])
            r_values.append(r)
            for key, val in zip(["SS", "SD", "DS", "DD"], parts[1:]):
                rdf[key].append(float(val) * N)

        print(f"Read RDF file '{filename}' with normalization count: {N}", file=sys.stderr)
        return N, np.array(r_values), {k: np.array(v) for k, v in rdf.items()}

    except Exception as e:
        print(f"Error reading RDF file '{filename}': {e}")
        return 0, [], {}

def classify_pair(h1, u1, h2, u2):
    same_hand = (h1 == h2)
    same_dir = (u1 == u2)
    if same_hand and same_dir:
        return "SS"
    elif same_hand and not same_dir:
        return "SD"
    elif not same_hand and same_dir:
        return "DS"
    else:
        return "DD"

def read_structure_file(filename, N, rVal, rdf):
	try:
		with open(filename, 'r') as f:
			lines = f.readlines()
	
		box_size_x, box_size_y = map(float, lines[0].strip().split())
		data_lines = lines[1:]
	
		x_vals, y_vals = [], []
		handedness, direction = [], []
	
		for line in data_lines:
			parts = line.strip().split()
			if len(parts) != 7:
				continue
			x = float(parts[0])
			y = float(parts[1])
			h = int(parts[4])
			d = int(parts[5])
			x_vals.append(x)
			y_vals.append(y)
			handedness.append(h)
			direction.append(d)
	
		x_vals = np.array(x_vals)
		y_vals = np.array(y_vals)
		handedness = np.array(handedness)
		direction = np.array(direction)
	
		nparts = len(x_vals)
		nbins = len(rVal)
	
		for i in range(nparts):
			for j in range(i + 1, nparts):
				dx = x_vals[i] - x_vals[j]
				dy = y_vals[i] - y_vals[j]
	
				# Apply minimum image convention
				dx -= box_size_x * round(dx / box_size_x)
				dy -= box_size_y * round(dy / box_size_y)
	
				r = np.sqrt(dx * dx + dy * dy)
	
				# Find bin
				for k in range(nbins - 1):
					if rVal[k] <= r < rVal[k + 1]:
						category = classify_pair(handedness[i], direction[i], handedness[j], direction[j])
						rdf[category][k] += 2  # count both (i,j) and (j,i)
						break
	
			N += 1  # one structure processed
		return N, rdf
	
	except Exception as e:
		print(f"Error processing structure file '{filename}': {e}")
		return N, rdf

def write_rdf(filename, N, rVal, rdf):
    try:
        with open(filename, 'w') as f:
            f.write(f"{N}\n")
            for i in range(len(rVal)):
                values = [rVal[i]] + [rdf[key][i] / max(N, 1) for key in ["SS", "SD", "DS", "DD"]]
                f.write("%12.8f %12.8E %12.8E %12.8E %12.8E\n" % tuple(values))
        print(f"Written RDF data to '{filename}'", file=sys.stderr)
    except Exception as e:
        print(f"Error writing RDF file '{filename}': {e}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Radial Distribution Function Tool (SS/SD/DS/DD, single normalization)")
    parser.add_argument("-init", nargs=3, metavar=('Rmin', 'Rmax', 'dR'), type=float,
                        help="Initialize RDF with given Rmin, Rmax, and dR")
    parser.add_argument("-read", nargs='?', const="rdf_bins.txt", metavar='filename',
                        help="Read and continue from existing RDF file")
    parser.add_argument("-file", metavar='structure_file',
                        help="Structure file with x/y, handedness, and direction info")
    parser.add_argument("-out", metavar='output_file', default='rdf_out.txt',
                        help="Output file for RDF data (default: rdf_out.txt)")

    args = parser.parse_args()

    N = 0
    rVal = []
    rdf = {}

    if args.init:
        rmin, rmax, dr = args.init
        N, rVal, rdf = init_mode(rmin, rmax, dr)
    elif args.read:
        N, rVal, rdf = read_mode(args.read)
    else:
        print("Error: You must specify either -init or -read mode.")
        sys.exit(1)

    if args.file:
        N, rdf = read_structure_file(args.file, N, rVal, rdf)

    write_rdf(args.out, N, rVal, rdf)

if __name__ == "__main__":
    main()
