import os, sys
import numpy as np
import pandas as pd

bobutils_dir = os.path.realpath("../code")
print(f"bobutils_dir={bobutils_dir}")
sys.path.append(bobutils_dir)

import bobutils.utils as bu
# Test harness function
def run_test_vectors(test_vectors, flux_cube, variance_cube):
    """
    Run test vectors to validate the fractional pixel mask on flux and variance cubes.
    
    :param test_vectors: A list of dictionaries with test vector data
    :param flux_cube: The 3D numpy array representing the flux
    :param variance_cube: The 3D numpy array representing the variance
    :return: A list with test results
    """
    test_results = []
    for test_vector in test_vectors:
        x_min, x_max, y_min, y_max = test_vector['x_min'], test_vector['x_max'], test_vector['y_min'], test_vector['y_max']
        expected_flux = test_vector['expected_flux']
        expected_var = test_vector['expected_var']
        
        # Create the mask
        mask = create_fractional_mask_corrected(x_min, x_max, y_min, y_max, flux_cube.shape[1:])
        
        # Apply the mask to the cubes
        calculated_flux, calculated_var = apply_fractional_mask_to_cubes(flux_cube, variance_cube, mask)
        
        # Compare the results
        flux_pass = np.allclose(calculated_flux, expected_flux, atol=1e-5)
        var_pass = np.allclose(calculated_var, expected_var, atol=1e-5)
        
        # Compile the results
        test_results.append({
            'Test Vector': test_vector,
            'Flux Pass': flux_pass,
            'Variance Pass': var_pass,
            'Calculated Flux': calculated_flux,
            'Calculated Variance': calculated_var
        })
        
    return test_results

# Example usage of the functions (you would replace this with actual data and test vectors)
# Create sample cubes
np.random.seed(0)  # Seed for reproducibility
flux_cube_checkerboard = np.indices((10, 10))[0] % 2
variance_cube_lower_diag = np.tril(np.ones((10, 10)))

# Extend the 2D patterns to 3D cubes
flux_cube_checkerboard = np.repeat(flux_cube_checkerboard[np.newaxis, :, :], 5, axis=0)
variance_cube_lower_diag = np.repeat(variance_cube_lower_diag[np.newaxis, :, :], 5, axis=0)

# Create a set of test vectors with expected outputs (this is an example, replace with actual test data)
test_vectors_with_expected_outputs = [
    {'x_min': 1.5, 'x_max': 3.5, 'y_min': 1.5, 'y_max': 3.5, 'expected_flux': np.array([2, 2, 2, 2, 2]), 'expected_var': np.array([4, 4, 4, 4, 4])},
    # ... more test vectors
]

# Run the test vectors
test_results = run_test_vectors(test_vectors_with_expected_outputs, flux_cube_checkerboard, variance_cube_lower_diag)

# Output the test results
for i, result in enumerate(test_results):
    print(f"Test Vector {i+1}: Flux {'PASS' if result['Flux Pass'] else 'FAIL'}, Variance {'PASS' if result['Variance Pass'] else 'FAIL'}")
