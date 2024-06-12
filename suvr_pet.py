import nibabel as nib
import numpy as np

def load_nifti_image(file_path):
    """Load a NIfTI image and return the image data and affine matrix."""
    img = nib.load(file_path)
    return img.get_fdata(), img.affine

def calculate_suvr(pet_data, mri_data, target_mask, reference_mask):
    """Calculate the SUVR given PET data, MRI data, target mask, and reference mask."""
    target_uptake = np.mean(pet_data[target_mask>0.5])
    reference_uptake = np.mean(pet_data[reference_mask>0.5])
    suvr = target_uptake / reference_uptake
    return suvr


# input files
path_to_pet_image = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri.nii.gz'
path_to_mri_image = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.nii.gz'
path_to_target_mask = '/home/ajoshi/Projects/PETBrainSuite/test_data/gm_mask.nii.gz'
path_to_reference_mask = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.cortex.dewisp.mask.nii.gz'

# output files
path_to_suvr = '/home/ajoshi/Projects/PETBrainSuite/test_data/suvr.nii.gz'



# Load PET and MRI images
pet_data, pet_affine = load_nifti_image(path_to_pet_image)
mri_data, mri_affine = load_nifti_image(path_to_mri_image)

# Define the target and reference masks (this can be done using MRI image-based segmentation)
# Here, assuming masks are already defined and loaded as NIfTI images
target_mask, _ = load_nifti_image(path_to_target_mask)
reference_mask, _ = load_nifti_image(path_to_reference_mask)

# Calculate SUVR
suvr = calculate_suvr(pet_data, mri_data, target_mask, reference_mask)
print("SUVR:", suvr)
