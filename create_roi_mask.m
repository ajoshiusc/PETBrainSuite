function create_roi_mask(labeled_image_path, roi_list, output_mask_path)
    % Read the labeled brain image
    labeled_image = niftiread(labeled_image_path);
    nii_info = niftiinfo(labeled_image_path);
    
    % Initialize the mask
    mask = zeros(size(labeled_image));
    
    % Iterate over the list of ROIs and create the mask
    for i = 1:length(roi_list)
        roi_value = roi_list(i);
        mask(labeled_image == roi_value) = 1;
    end
    
    % Save the mask as a NIfTI file
    nii_info.Datatype='uint8';
    niftiwrite(uint8(mask), output_mask_path, nii_info);
    
    disp(['Mask file saved to ', output_mask_path]);
end
