import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec


bobutils_dir = os.path.realpath("../code/bobutils")
print(f"bobutils_dir={bobutils_dir}")
sys.path.append(bobutils_dir)

import utils as bu

# Define a function to generate a rotated rectangle around the center point
def generate_rotated_rectangle(center, width, height, angle):
    """
    Generate the corners of a rotated rectangle.

    :param center: The (x, y) coordinates of the center of the rectangle
    :param width: The width of the rectangle
    :param height: The height of the rectangle
    :param angle: The rotation angle in degrees
    :return: An array of the corner points of the rotated rectangle
    """
    # Define the unrotated corners of the rectangle
    c_x, c_y = center
    corners = np.array([
        [c_x - width / 2, c_y - height / 2],
        [c_x + width / 2, c_y - height / 2],
        [c_x + width / 2, c_y + height / 2],
        [c_x - width / 2, c_y + height / 2]
    ])
    print(f"corners={corners}")
    # Rotate the corners around the center point
    theta = np.radians(angle)  # Convert angle to radians
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    rotated_corners = (rotation_matrix @ (corners - center).T).T + center
    print(f"rotated_corners={rotated_corners}")
    return rotated_corners

def generate_rotated_rectangle2(center, width, height, angle):
    """
    Generate the corners of a rotated rectangle.
    
    :param center: The (x, y) coordinates of the rectangle's center
    :param width: The width of the rectangle
    :param height: The height of the rectangle
    :param angle: The rotation angle in degrees, counter-clockwise from the x-axis
    :return: A numpy array of shape (4, 2) containing the x and y coordinates of the rectangle's corners
    """
    # Define the unrotated rectangle corners relative to the origin
    unrotated_corners = np.array([
        [-width / 2, -height / 2],
        [width / 2, -height / 2],
        [width / 2, height / 2],
        [-width / 2, height / 2]
    ])
    
    # Rotation matrix
    theta = np.radians(angle)
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    
    # Rotate the corners and translate them to the center point
    rotated_corners = np.dot(unrotated_corners, rotation_matrix) + center
    
    return rotated_corners

# Update the plotting function to use a colormap that accurately reflects the weight values
def plot_fractional_extraction_box_with_mask(mask, box_corners, ax, title=None):
    """ 
    we assume that the mask is the same size as the x,y spatial extent of the cube 
    so instead of doing a flux_cut and var_cut like we've been doing, we'll just apply the mask directly to the cube to come up with the flux and variance spectra.
    """
    # Define a colormap from white to black
    cmap = plt.cm.gray_r

    # Normalize the colormap to the range of weights in the mask
    norm = mcolors.Normalize(vmin=0, vmax=1)

    # Plot the pixel grid and weights
    for i in range(mask.shape[1]):
        for j in range(mask.shape[0]):
            # Only label and color non-zero weights
            if mask[j, i] != 0:
                weight = mask[j, i]
                color = cmap(norm(weight))
                text_color = 'white' if weight > 0.5 else 'black'
                ax.text(i + 0.5, j + 0.5, f"{weight:.2f}", ha='center', va='center', color=text_color)
                ax.add_patch(plt.Rectangle((i, j), 1, 1, facecolor=color))

    # Plot the extraction box outline in red
    ax.add_patch(plt.Polygon(box_corners, closed=True, edgecolor='r', fill=False))

    # Plot the dotted lines for the bounding box edges
    for point in box_corners:
        ax.axhline(y=point[1], color='r', linestyle=':', linewidth=1)
        ax.axvline(x=point[0], color='r', linestyle=':', linewidth=1)
        ax.text(point[0], -0.5, f"{point[0]:.2f}", ha='center', va='center', color='r')
        ax.text(-0.5, point[1], f"{point[1]:.2f}", ha='center', va='center', color='r', rotation=90)

    # Set the plot limits
    xlims = [0,10]
    ylims = [0,10]
    ax.set_xlim(0, mask.shape[1])
    ax.set_ylim(0, mask.shape[0])

    # Add grid, ticks, and set aspect ratio
    ax.grid(True, which='both', color='k', linestyle='-', linewidth=1)
    ax.set_xticks(range(mask.shape[1] + 1))
    ax.set_yticks(range(mask.shape[0] + 1))
    ax.set_aspect('equal')
    ax.set_axisbelow(True)
    ax.set_facecolor('white')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if title is not None: ax.set_title(title, fontsize=14, pad=20, fontweight='bold')


    return ax

def calculate_bounding_box(mask):
    """
    Calculate the bounding box for non-zero weighted pixels in a mask.

    :param mask: A 2D numpy array representing the mask with fractional weights
    :return: A tuple containing the bounding box coordinates (x_min, x_max, y_min, y_max)
    """
    # Find indices of non-zero elements
    non_zero_indices = np.nonzero(mask)
    y_indices, x_indices = non_zero_indices

    # Determine the bounding box coordinates
    x_min, x_max = x_indices.min(), x_indices.max()
    y_min, y_max = y_indices.min(), y_indices.max()

    # Return the bounding box coordinates
    # Adding 1 to x_max and y_max to ensure the max index is included
    return x_min, x_max + 1, y_min, y_max + 1

def find_minimal_mask_shape(rotated_corners):
    """
    Find the minimal mask shape that contains the rotated rectangle, adjusting for cases where
    the corners align with the pixel grid.
    
    :param rotated_corners: The corners of the rotated rectangle
    :return: The minimal mask shape as a tuple (height, width)
    """
    # Find the axis-aligned bounding box of the rotated rectangle
    x_min, y_min = np.min(rotated_corners, axis=0)
    x_max, y_max = np.max(rotated_corners, axis=0)
    
    # Calculate the shape as the ceiling of the differences
    # Add 1 to ensure the rectangle fits within the mask if the max coordinates are not integers
    height = int(np.ceil(y_max)) - int(np.floor(y_min))
    width = int(np.ceil(x_max)) - int(np.floor(x_min))
    
    return (height, width)

from shapely.geometry import Polygon

def calc_overlap_area(box_corners, pixel_corners):
    """
    Calculate the area of overlap between the extraction box and a pixel.
    
    :param box_corners: Corners of the extraction box
    :param pixel_corners: Corners of the pixel
    :return: The overlap area
    """
    # Create shapely polygons for the box and the pixel
    box_poly = Polygon(box_corners)
    pixel_poly = Polygon(pixel_corners)
    
    # Calculate the area of overlap
    overlap_area = box_poly.intersection(pixel_poly).area
    return overlap_area


# Let's redefine the create_fractional_mask_rotated function to correctly calculate the weights
def create_fractional_mask_rotated(box_corners, shape):
    """
    Create a fractional mask for a given bounding box with subpixel accuracy, handling rotation.
    
    :param box_corners: The corners of the bounding box (rotated)
    :param shape: The shape of the 2D array to apply the mask to
    :return: A 2D numpy array representing the mask
    """
    mask = np.zeros(shape)
    
    # Iterate over the coordinates of the pixels
    for y in range(shape[0]):
        for x in range(shape[1]):
            # Define the corners of the current pixel
            pixel_corners = np.array([
                [x, y],
                [x+1, y],
                [x+1, y+1],
                [x, y+1]
            ])
            # Calculate the area of overlap and set it as the weight for the pixel
            mask[y, x] = calc_overlap_area(box_corners, pixel_corners)
    
    return mask


def run_test(ax, title, rotation_degrees):
    # Use the function to generate a rectangle rotated by 15 degrees
    center_point = [5.2, 4.87]
    rect_width = np.sqrt(3)  # Width such that the area is greater than 3 square units
    rect_height = 3  # Height arbitrarily chosen
    mask_shape = (10,10)
    rotated_rect_corners = generate_rotated_rectangle(center_point, rect_width, rect_height, rotation_degrees)
    # Create the mask for the rotated rectangle
    rotated_rect_mask = create_fractional_mask_rotated(rotated_rect_corners, mask_shape)
    ax = plot_fractional_extraction_box_with_mask(rotated_rect_mask, rotated_rect_corners, ax, title)


def main():
    # Create the plot with the accurate colormap for the rotated rectangle
    fig = plt.figure(figsize=(10,10), layout="constrained", dpi=150)
    gs = gridspec.GridSpec(1, 2, figure=fig)
    ax_plain = fig.add_subplot(gs[0, 0])
    ax_rotated = fig.add_subplot(gs[0, 1])

    run_test(ax_plain, "Fractional Mask for a Rectangular Extraction Box", 0)
    run_test(ax_rotated, "Fractional Mask for a Rotated Extraction Box", 15)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()