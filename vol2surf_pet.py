import nibabel as nib
import numpy as np
import scipy.ndimage as ndimage
from nibabel.freesurfer import read_geometry
from dfsio import readdfs
from surfproc import view_patch_vtk, patch_color_attrib, smooth_patch, view_patch, view_patch_vtk


def load_nifti_image(file_path):
    """Load a NIfTI image and return the image data and affine matrix."""
    img = nib.load(file_path)
    return img.get_fdata(), img.affine

def load_cortical_surface(surface_path):
    """Load a cortical surface file and return the vertex coordinates and faces."""
    s = readdfs(surface_path)
    coords = s.vertices.copy()
    faces = s.faces.copy()


    return coords, faces

def sample_pet_to_surface(pet_data, pet_affine, surface_coords):
    """Sample PET data onto the cortical surface."""
    # Transform surface coordinates to voxel space
    surface_coords_voxel = surface_coords #nib.affines.apply_affine(np.linalg.inv(pet_affine), surface_coords)

    # Sample PET values at surface vertices
    sampled_values = ndimage.map_coordinates(pet_data, surface_coords_voxel.T, order=1)
    return sampled_values

def smooth_surface_data(data, faces, iterations=10):
    """Smooth data on the cortical surface using iterative averaging."""
    smoothed_data = data.copy()
    for _ in range(iterations):
        for i, face in enumerate(faces):
            smoothed_data[face] = smoothed_data[face].mean()
    return smoothed_data

# File paths
pet_image_path = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet_pvc.nii.gz'
cortical_surface_path = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.left.mid.cortex.svreg.dfs'

# Load PVC-corrected PET volume
pet_data, pet_affine = load_nifti_image(pet_image_path)

# Load cortical surface
surface_coords, surface_faces = load_cortical_surface(cortical_surface_path)

# Sample PET data onto the cortical surface
sampled_pet_values = sample_pet_to_surface(pet_data, pet_affine, surface_coords)

# Smooth the sampled PET data
smoothed_pet_values = smooth_surface_data(sampled_pet_values, surface_faces, iterations=10)

# visualize the smoothed PET values on the cortical surface and save to a png file using nilearn
import nilearn.plotting as nlp

#nlp.view_surf(surface_coords, surface_faces, smoothed_pet_values, cmap='viridis') #, symmetric_cmap=False, colorbar=True, threshold=None, bg_map=None, vmax=None, figure=None, title='Smoothed PET values on cortical surface', output_file='smoothed_pet_surface.png', display_mode='auto', view='lateral')

print("Smoothed PET values visualized on the cortical surface and saved to 'smoothed_pet_surface.png'")

# Alternatively, visualize the smoothed PET values on the cortical surface using matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

""" fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(surface_coords[:, 0], surface_coords[:, 1], surface_coords[:, 2], triangles=surface_faces, cmap='viridis', alpha=0.5)
sc = ax.scatter(surface_coords[:, 0], surface_coords[:, 1], surface_coords[:, 2], c=smoothed_pet_values, cmap='viridis')
plt.colorbar(sc)
plt.savefig('smoothed_pet_surface.png')
plt.show() """


class surface:
    def __init__(self, coords, faces):
        self.vertices = coords
        self.faces = faces

s = surface(surface_coords, surface_faces)
s.attributes = smoothed_pet_values #sampled_pet_values

s = patch_color_attrib(s, s.attributes, cmap='jet', clim=[0, 1000])

view_patch_vtk(s)



# Save smoothed PET values to a file (e.g., text file)
np.savetxt('smoothed_pet_values.txt', smoothed_pet_values)

print("Smoothed PET values saved to 'smoothed_pet_values.txt'")
