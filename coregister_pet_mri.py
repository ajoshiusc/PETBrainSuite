import SimpleITK as sitk

def load_nifti_image(file_path):
    """Load a NIfTI image and return the SimpleITK image object."""
    return sitk.ReadImage(file_path)

def save_nifti_image(sitk_image, output_path):
    """Save a SimpleITK image object as a NIfTI file."""
    sitk.WriteImage(sitk_image, output_path)

def coregister_images(fixed_image_path, moving_image_path, output_path):
    """Coregister the moving image to the fixed image."""
    # Load fixed (MRI) and moving (PET) images
    fixed_image = load_nifti_image(fixed_image_path)
    moving_image = load_nifti_image(moving_image_path)

    # Initial alignment of the two volumes
    initial_transform = sitk.CenteredTransformInitializer(fixed_image, 
                                                          moving_image, 
                                                          sitk.Euler3DTransform(), 
                                                          sitk.CenteredTransformInitializerFilter.GEOMETRY)

    # Set up the registration framework
    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, 
                                                      numberOfIterations=100, 
                                                      convergenceMinimumValue=1e-6, 
                                                      convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Setup for the multi-resolution framework.
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Don't optimize in-place, we would possibly like to run this cell multiple times.
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    # Perform the registration.
    final_transform = registration_method.Execute(sitk.Cast(fixed_image, sitk.sitkFloat32), 
                                                  sitk.Cast(moving_image, sitk.sitkFloat32))

    # Resample the moving image to the fixed image space.
    moving_resampled = sitk.Resample(moving_image, fixed_image, final_transform, sitk.sitkLinear, 0.0, moving_image.GetPixelID())

    # Save the registered image.
    save_nifti_image(moving_resampled, output_path)
    print(f"Coregistered image saved to {output_path}")

# File paths
fixed_image_path = 'path_to_mri_image.nii'
moving_image_path = 'path_to_pet_image.nii'
output_path = 'path_to_registered_pet_image.nii'

# Coregister PET to MRI
coregister_images(fixed_image_path, moving_image_path, output_path)
