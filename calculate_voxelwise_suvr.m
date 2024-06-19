function calculate_voxelwise_suvr(pet_image_path, reference_mask_path, output_suvr_path)
    % Read the PET image
    pet_data = niftiread(pet_image_path);
    nii_info = niftiinfo(pet_image_path);
    
    % Read the reference region mask
    reference_mask = niftiread(reference_mask_path);
    
    % Calculate the mean uptake in the reference region
    reference_mean = mean(pet_data(reference_mask > 0));
    
    % Calculate voxelwise SUVR
    suvr_data = pet_data / reference_mean;
    
    % Save the SUVR image as a NIfTI file
    niftiwrite(suvr_data, output_suvr_path, nii_info);
    
    disp(['Voxelwise SUVR image saved to ', output_suvr_path]);
end
