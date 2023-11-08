# code and test vectors
import numpy as np
import pandas as pd

# Function to create the fractional mask
def create_fractional_mask_corrected(x_min, x_max, y_min, y_max, shape):
    """
    Create a fractional mask for a given bounding box with subpixel accuracy.
    
    :param x_min: The minimum x-coordinate of the bounding box
    :param x_max: The maximum x-coordinate of the bounding box
    :param y_min: The minimum y-coordinate of the bounding box
    :param y_max: The maximum y-coordinate of the bounding box
    :param shape: The shape of the 2D array to apply the mask to
    :return: A 2D numpy array representing the mask
    """
    y, x = np.ogrid[:shape[0], :shape[1]]
    mask = np.zeros(shape)
    mask[y, x] = ((x >= x_min) & (x < x_max) & (y >= y_min) & (y < y_max)).astype(float)
    
    # Add fractional parts of the pixels
    mask[y, x] += (x == x_min) * (1 - (x_min % 1)) * (y >= y_min) * (y < y_max)
    mask[y, x] += (x == x_max) * (x_max % 1) * (y >= y_min) * (y < y_max)
    mask[y, x] += (y == y_min) * (1 - (y_min % 1)) * (x >= x_min) * (x < x_max)
    mask[y, x] += (y == y_max) * (y_max % 1) * (x >= x_min) * (x < x_max)
    
    # Correct the corners
    mask[y == y_min, x == x_min] = (1 - (x_min % 1)) * (1 - (y_min % 1))
    mask[y == y_min, x == x_max] = (x_max % 1) * (1 - (y_min % 1))
    mask[y == y_max, x == x_min] = (1 - (x_min % 1)) * (y_max % 1)
    mask[y == y_max, x == x_max] = (x_max % 1) * (y_max % 1)
    
    return mask

# Function to apply the mask to the flux and variance cubes
def apply_fractional_mask_to_cubes(flux_cube, variance_cube, mask):
    """
    Apply the fractional pixel mask to the flux and variance cubes.
    
    :param flux_cube: The 3D numpy array representing the flux
    :param variance_cube: The 3D numpy array representing the variance
    :param mask: The 2D numpy array representing the fractional mask
    :return: Tuple of the masked flux and variance arrays
    """
    # Apply the mask to each layer of the flux and variance cubes
    masked_flux = flux_cube * mask[np.newaxis, :, :]
    masked_variance = variance_cube * mask[np.newaxis, :, :]
    
    # Sum the masked arrays across the spatial dimensions to get the spectra
    flux_spectrum = masked_flux.sum(axis=(1, 2))
    variance_spectrum = masked_variance.sum(axis=(1, 2))
    
    return flux_spectrum, variance_spectrum

