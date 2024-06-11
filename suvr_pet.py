import nibabel as nib
import numpy as np

def load_nifti_image(file_path):
    """Load a NIfTI image and return the image data and affine matrix."""
    img = nib.load(file_path)
    return img.get_fdata(), img.affine

def calculate_suvr(pet_data, mri_data, target_mask, reference_mask):
    """Calculate the SUVR given PET data, MRI data, target mask, and reference mask."""
    target_uptake = np.mean(pet_data[target_mask])
    reference_uptake = np.mean(pet_data[reference_mask])
    suvr = target_uptake / reference_uptake
    return suvr

# Load PET and MRI images
pet_data, pet_affine = load_nifti_image('path_to_pet_image.nii')
mri_data, mri_affine = load_nifti_image('path_to_mri_image.nii')

# Define the target and reference masks (this can be done using MRI image-based segmentation)
# Here, assuming masks are already defined and loaded as NIfTI images
target_mask, _ = load_nifti_image('path_to_target_mask.nii')
reference_mask, _ = load_nifti_image('path_to_reference_mask.nii')

# Calculate SUVR
suvr = calculate_suvr(pet_data, mri_data, target_mask, reference_mask)
print("SUVR:", suvr)
