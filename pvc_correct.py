import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter

def load_nifti_image(file_path):
    """Load a NIfTI image and return the image data and affine matrix."""
    img = nib.load(file_path)
    return img.get_fdata(), img.affine

def muller_gaertner_pvc(pet_data, mri_data, gray_matter_mask, white_matter_mask, fwhm):
    """Perform Müller-Gärtner PVC on PET data using MRI segmentation masks."""
    # Smoothing the PET image to match spatial resolution
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    smoothed_pet = gaussian_filter(pet_data, sigma=sigma)

    # Calculate mean uptake in white matter (assumed to be constant across the brain)
    mean_white_matter_uptake = np.mean(smoothed_pet[white_matter_mask>0.5])

    # Calculate corrected PET uptake for gray matter
    corrected_pet = np.zeros(pet_data.shape)
    corrected_pet[gray_matter_mask>0.5] = (pet_data[gray_matter_mask>0.5] - mean_white_matter_uptake) / (1.0001 - gray_matter_mask.astype(float)[gray_matter_mask>0.5])

    return corrected_pet


# input files
path_to_pet_image = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri.nii.gz'
path_to_mri_image = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.nii.gz'
path_to_white_matter_mask = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.cortex.dewisp.mask.nii.gz'

path_to_gm_frac = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.GM_frac.nii.gz'


# output files
path_to_gray_matter_mask = '/home/ajoshi/Projects/PETBrainSuite/test_data/gm_mask.nii.gz'
path_to_corrected_pet_image = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet_pvc.nii.gz'


# read the GM fraction, threshold it at 0.5, and write it to a mask
gm_frac = nib.load(path_to_gm_frac).get_fdata()
gm_mask = gm_frac > 0.5
gm_mask_img = nib.Nifti1Image(gm_mask.astype(np.uint8), nib.load(path_to_mri_image).affine)
nib.save(gm_mask_img, path_to_gray_matter_mask)
print("Gray matter mask saved.")




# Load PET and MRI images
pet_data, pet_affine = load_nifti_image(path_to_pet_image)
mri_data, mri_affine = load_nifti_image(path_to_mri_image)

# Load segmentation masks (gray matter and white matter)
gray_matter_mask, _ = load_nifti_image(path_to_gray_matter_mask)
white_matter_mask, _ = load_nifti_image(path_to_white_matter_mask)

# Perform Müller-Gärtner PVC
fwhm = 6.0  # Example Full Width at Half Maximum (FWHM) in mm, adjust based on your PET scanner
corrected_pet = muller_gaertner_pvc(pet_data, mri_data, gray_matter_mask, white_matter_mask, fwhm)

# Save corrected PET image
corrected_pet_img = nib.Nifti1Image(corrected_pet, pet_affine)
nib.save(corrected_pet_img, path_to_corrected_pet_image)
print("PVC corrected PET image saved.")
