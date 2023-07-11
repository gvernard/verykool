import os
import sys
import time
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm
from pprint import pprint

from coolest.api.plotting import ModelPlotter, MultiModelPlotter
from coolest.api import util


def get_coolest_template(code_choice, model_choice):
    json_files = glob.glob(os.path.join('coolest_results', code_choice, model_choice, '*.json'))
    file_no_ext = os.path.splitext(json_files[0])[0]
    if file_no_ext[-6:] == '_pyAPI':
        file_no_ext = file_no_ext[:-6]
    return file_no_ext

def get_source_title(coolest, galaxy_idx=2):
    source_model = coolest.lensing_entities[galaxy_idx].light_model
    profile_names = [profile.type for profile in source_model]
    return str(profile_names)





vkl_path = sys.argv[1] #"/home/giorgos/myData/VKL_runs/sourcemag22_unfiltered/test/VKL_coolest/coolest_vkl"
coolest_vkl = util.get_coolest_object(vkl_path+"VKL_coolest/coolest_vkl")


coordinates = util.get_coordinates(coolest_vkl)
x, y = coordinates.pixel_coordinates
print(coordinates.plt_extent)
coordinates_src = coordinates.create_new_coordinates(pixel_scale_factor=0.1,grid_shape=(1.8, 1.8))
x_src, y_src = coordinates_src.pixel_coordinates
print(coordinates_src.plt_extent)


coolest_vkl.lensing_entities[-1].light_model


# initialize the plotter
plotter = MultiModelPlotter([coolest_vkl],
                            coolest_directories=[os.path.dirname(vkl_path)]
                           )
plotter_irreg = ModelPlotter(coolest_vkl, coolest_directory=os.path.dirname(vkl_path))

norm = None # Normalize(-1e-4)
neg_values = False


# create the figure
fig, axes = plt.subplots(5, 1, figsize=(3, 12))

# multi-plotter panels
plotter.plot_surface_brightness(
    [axes[0]], 
    titles=[
        get_source_title(coolest_vkl), 
    ],
    kwargs_light=dict(entity_selection=[[2]]),
    norm=norm,
    neg_values_as_bad=neg_values,
    coordinates=coordinates_src)

# single plotter panels
plotter_irreg.plot_surface_brightness(
    axes[1], title=get_source_title(coolest_vkl),
    kwargs_light=dict(entity_selection=[2]),
    norm=norm, 
    neg_values_as_bad=neg_values,
    extent=coordinates_src.extent)

plotter.plot_model_image(
    [axes[2]],
    supersampling=5, convolved=True,
    kwargs_source=dict(entity_selection=[[2]]),
    kwargs_lens_mass=dict(entity_selection=[[0, 1]]),
    norm=norm)
#axes[2].set_xlim(-2.5, 2.5);axes[2].set_ylim(-2.5, 2.5)
#axes[2, 1].set_xlim(-2.5, 2.5);axes[2, 1].set_ylim(-2.5, 2.5)
plotter.plot_model_residuals(
    [axes[3]],
    supersampling=5,
    kwargs_source=dict(entity_selection=[[2]]),
    kwargs_lens_mass=dict(entity_selection=[[0, 1]]),
    norm=norm)
plotter.plot_magnification(
    [axes[4]],
    kwargs_lens_mass=dict(entity_selection=[[0, 1]]),
    norm=norm)

#axes[2, 2].axis('off')
#axes[3, 2].axis('off')
#axes[4, 2].axis('off')

fig.tight_layout()
plt.savefig("coolest_plot.png")
