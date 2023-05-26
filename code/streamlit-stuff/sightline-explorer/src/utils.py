import numpy as np
import pandas as pd
from PIL import Image, ImageEnhance

def show_image_stats(label, wl_image):
    print(f"---- {label} ------")
    # Number of elements in the image
    num_elements = np.prod(wl_image.shape)

    # Number of NaN values
    num_nan = np.sum(np.isnan(wl_image))

    # Number of negative values
    num_negative = np.sum(wl_image < 0)

    # Minimum and maximum values
    min_val = np.nanmin(wl_image)
    max_val = np.nanmax(wl_image)

    # Mean and standard deviation
    mean_val = np.nanmean(wl_image)
    std_val = np.nanstd(wl_image)

    print("---- BEFORE CORRECTIONS ------")
    print(f"Number of elements: {num_elements}")
    print(f"Number of NaN values: {num_nan}")
    print(f"Number of negative values: {num_negative}")
    print(f"Min value: {min_val}")
    print(f"Max value: {max_val}")
    print(f"Mean value: {mean_val}")
    print(f"Standard deviation: {std_val}")

def make_image_corrections(wl_image, contrast, brightness, sharpness, scale_factor):
    # Replace NaN values with 0
    wl_image = np.nan_to_num(wl_image)

    # Scale the data to the range 0-1
    wl_image = (wl_image - np.min(wl_image)) / (np.max(wl_image) - np.min(wl_image))

    # Flip the y-axis
    wl_image = np.flipud(wl_image)

    # Invert and convert to PIL image
    wl_image = Image.fromarray((wl_image * 255).astype(np.uint8))

    # original size
    original_size = wl_image.size

    # new size
    new_size = [dimension * scale_factor for dimension in original_size]

    # resize the image
    wl_image = wl_image.resize(new_size, Image.NEAREST)

    # Adjust contrast
    enhancer = ImageEnhance.Contrast(wl_image)
    wl_image = enhancer.enhance(contrast)

    # Adjust brightness
    enhancer = ImageEnhance.Brightness(wl_image)
    wl_image = enhancer.enhance(brightness)   

    # Adjust brightness
    enhancer = ImageEnhance.Sharpness(wl_image)
    wl_image = enhancer.enhance(sharpness)  

    return wl_image

# handline dataframes

# add/append a row to the specified dataframe 
def append_row(df, new_row):
    new_row_df = pd.DataFrame([new_row])
    df = pd.concat([df, new_row_df], ignore_index=True)
    return df

# remove a row from the specified dataframe using the provided index
def remove_row_with_index(df, index):
    df = df.drop(index)
    return df

# remove a row (or rows) from the specified dataframe using a column label and a key value
def remove_row_with_key(df, column, key):
    df = df[df[column] != key]
    return df