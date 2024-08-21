function create_roi_mask(labeled_image_path, roi_list,pvc_label_path, output_mask_path)
    % Read the labeled brain image
    labeled_image = niftiread(labeled_image_path);
    nii_info = niftiinfo(labeled_image_path);    
    % Initialize the mask
    mask = zeros(size(labeled_image));

    pvc_label= niftiread(pvc_label_path);

    % pvc_label == 2 is Gray Matter
    % Iterate over the list of ROIs and create the mask
    for i = 1:length(roi_list)
        roi_value = roi_list(i);
        mask(labeled_image == roi_value & pvc_label == 2) = 1;
    end
     
    % Save the mask as a NIfTI file
    nii_info.Datatype='uint8';
    niftiwrite(uint8(mask), output_mask_path, nii_info);
    
    disp(['Mask file saved to ', output_mask_path]);
end
