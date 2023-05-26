import streamlit as st
import streamlit.components.v1 as components
from streamlit_drawable_canvas import st_canvas

import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import os

from kcwitools import io as kcwi_io
from kcwitools import spec as kcwi_s
from kcwitools import utils as kcwi_u
from kcwitools.image import build_whitelight
import utils as utils

matplotlib.use('Agg')

st.set_page_config(
    page_title="sightline-explorer",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="collapsed",
    menu_items={
        'Get Help': 'https://www.extremelycoolapp.com/help',
        'Report a bug': "https://www.extremelycoolapp.com/bug",
        'About': "# This is a header. This is an *extremely* cool app!"
    }
)
# Create file uploaders for the flux and variance data cubes
# flux_file = st.file_uploader("Upload flux data cube", type=['fits', 'fits.gz'])
# variance_file = st.file_uploader("Upload variance data cube", type=['fits', 'fits.gz'])
# to use this we would need to implement the streaming version of load_file found in kcwi_tools.
# so we go basic
# if flux_filename is None or var_filename is None:
#     st.error('Please upload both the flux and variance files.')
# if hdr is None or flux is None or var is None:
#     st.error('Unable to load data from the uploaded files.')

st.cache_data
def load_fits_files(flux_filename,var_filename):
    # load flux and variance fits file
    base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429/"
    # flux_filename = base_path+"/J1429_rb_flux.fits"
    # var_filename = base_path+"/J1429_rb_var.fits"

    hdr, flux = kcwi_io.open_kcwi_cube(base_path+flux_filename)
    wave = kcwi_u.build_wave(hdr)
    _, var = kcwi_io.open_kcwi_cube(base_path+var_filename)
    return hdr, wave, flux, var

hdr, wave, flux, var = load_fits_files("J1429_rb_flux.fits","J1429_rb_var.fits")


# print(wl_image)
wl_image=build_whitelight(hdr, flux, minwave=4658, maxwave=4665)
# utils.show_image_stats("BEFORE CORRECTIONS", wl_image)
wl_image, wl_image_pil = utils.make_image_corrections(wl_image)
# utils.show_image_stats("AFTER CORRECTIONS", wl_image)

# Create a figure and axes
fig, ax = plt.subplots()
ax.imshow(wl_image,interpolation="nearest",cmap="gray",vmin=0)

# how about a little layout
image_area, sightlines_area, spectrum_area = st.columns([2, 1, 2])

with image_area:
   st.write("Flux")
   st.pyplot(fig)

with sightlines_area:
   st.write("Sightlines")
   st.image("https://static.streamlit.io/examples/dog.jpg")

with spectrum_area:
   st.write("Spectra")
   st.image("https://static.streamlit.io/examples/owl.jpg")



# Parameters for the interactive canvas
stroke_width = 3
stroke_color = "#ffffff"
bg_color = "#000000"


drawing_mode = st.sidebar.selectbox(
    "Drawing tool:", ("point", "freedraw", "line", "rect", "circle", "transform")
)
stroke_width = st.sidebar.slider("Stroke width: ", 1, 25, 3)
# if drawing_mode == 'point':
#     point_display_radius = st.sidebar.slider("Point display radius: ", 1, 25, 3)
stroke_color = st.sidebar.color_picker("Stroke color hex: ")
bg_color = st.sidebar.color_picker("Background color hex: ", "#eee")
bg_image = st.sidebar.file_uploader("Background image:", type=["png", "jpg"])
realtime_update = st.sidebar.checkbox("Update in realtime", True)



# Create the interactive canvas
canvas_result = st_canvas(
    fill_color="rgba(255, 165, 0, 0.3)",  # Fixed fill color with some opacity
    stroke_width=stroke_width,
    stroke_color=stroke_color,
    background_color=bg_color,
    background_image=wl_image_pil,
    update_streamlit=realtime_update,
    # height=200,
    # width=200,
    drawing_mode=drawing_mode,
    key="canvas",
)

# Display the coordinates of the selected point
if canvas_result.json_data is not None:
    st.write("Objects created:")
    st.json(canvas_result.json_data)

    # for path in canvas_result.json_data["objects"]:
    #     st.write(f"x: {path['left']}, y: {path['top']}")
