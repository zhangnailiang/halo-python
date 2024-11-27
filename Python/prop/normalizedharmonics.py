def normalizedharmonics(filepath, degree):
    """
    Reads the gravity model file and extracts the normalized spherical harmonic coefficients
    up to the specified degree.

    Parameters:
    - filepath: str, path to the gravity model file
    - degree: int, maximum degree of the harmonics

    Returns:
    - Cnm: 2D NumPy array, normalized cosine coefficients
    - Snm: 2D NumPy array, normalized sine coefficients
    """
    import numpy as np

    # Initialize coefficient arrays
    Cnm = np.zeros((degree + 1, degree + 1))
    Snm = np.zeros((degree + 1, degree + 1))

    with open(filepath, 'r') as file:
        for line in file:
            # Skip comments or empty lines
            line = line.strip()
            if not line or line.startswith(('#', '//', '%')):
                continue

            # Split the line into components
            parts = line.split()

            # Check if the line has enough parts
            if len(parts) < 4:
                continue

            try:
                # Parse degree and order
                n = int(parts[0])
                m = int(parts[1])

                if n > degree or m > degree:
                    continue  # Skip if beyond maximum degree or order

                # Parse coefficients
                Cnm_val = float(parts[2])
                Snm_val = float(parts[3])

                # Assign coefficients
                Cnm[n, m] = Cnm_val
                Snm[n, m] = Snm_val

            except ValueError:
                # Skip lines with invalid data
                continue

    return Cnm, Snm