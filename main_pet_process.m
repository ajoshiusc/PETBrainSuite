% Read the tissue fraction NIfTI files
clc;clear all;close all;restoredefaultpath;
addpath('/home/ajoshi/Projects/svreg/src');
addpath('/home/ajoshi/Projects/svreg/3rdParty');

% PET image coregistered to MRI
pet2mri_file = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri2_itksnap.nii.gz'; 

% PSF of PET image (from manufacturer)
xres = '6.0'; yres = '6.0'; zres = '6.0';

% MRI image should be processed using BrainSuite upto SVReg
mri_basefile = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri';
PETPVC_EXE = '/home/ajoshi/Projects/PETBrainSuite/PETPVC-1.2.10/bin/petpvc';

REF_ROIS = [900]% List of ROIS from atlas used as reference




% Output filenames
pvc_mask_file = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri2_itksnap_mask.nii.gz';
pet2mri_pvc_file = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri2_itksnap_pvc.nii.gz'; 
pet2surf_file = '/home/ajoshi/Projects/PETBrainSuite/test_data/pet2mri2_itksnap_pvc.mat';
pet_suvr_file = [mri_basefile,'.pet.suvr.nii.gz'];
pet_ref_mask_file = [mri_basefile,'.pet.ref.mask.nii.gz'];




tissueFracFile = [mri_basefile,'.pvc.frac.nii.gz'];



tissueFracData = niftiread(tissueFracFile);

grayMatterData = (3-tissueFracData).*(tissueFracData>=2) +(tissueFracData-1).*(tissueFracData<2).*(tissueFracData>=1);
whiteMatterData = (tissueFracData-2).*(tissueFracData>2);
CSFData = (2-tissueFracData).*(tissueFracData<=2).*(tissueFracData>1) + (tissueFracData).*(tissueFracData<=1);
CSFData(tissueFracData<=0)=1;
% Check if the dimensions match
if ~isequal(size(grayMatterData), size(whiteMatterData))
    error('Gray matter and white matter data must have the same dimensions.');
end

% Create a 4D array with the 4th dimension as channels
outputData = cat(4, grayMatterData, whiteMatterData,CSFData);

% Load the NIfTI information from one of the input files
niftiInfo = niftiinfo(tissueFracFile);

% Update the information for the 4D file
niftiInfo.ImageSize = size(outputData);
niftiInfo.PixelDimensions = [niftiInfo.PixelDimensions, 2]; % Add the 4th dimension

% Save the new 4D NIfTI file
niftiwrite(outputData, pvc_mask_file, niftiInfo);


cmd = [PETPVC_EXE,' -i ',pet2mri_file,' -m ',pvc_mask_file,' -o ',pet2mri_pvc_file,' --pvc IY -x ',xres,' -y ',yres,' -z ',zres];
system(cmd);



% Create mask file for reference rois
create_roi_mask([mri_basefile,'.svreg.label.nii.gz'], REF_ROIS, [mri_basefile,'.pet.ref.mask.nii.gz'])

% Calculate SUVR
calculate_voxelwise_suvr(pet2mri_pvc_file, pet_ref_mask_file,  pet_suvr_file)




%resample suvr to surface
resample2surf(mri_basefile,pet_suvr_file,pet2surf_file, 0)

a = load(pet2surf_file);

left_surf = readdfs([mri_basefile, '.left.mid.cortex.svreg.dfs']);
right_surf = readdfs([mri_basefile, '.right.mid.cortex.svreg.dfs']);

figure;
patch('faces',left_surf.faces,'Vertices',left_surf.vertices,'facevertexcdata',a.datal,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(-90,0);clim([0,3]);
camlight;

figure;
patch('faces',right_surf.faces,'Vertices',right_surf.vertices,'facevertexcdata',a.datar,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(90,0);clim([0,3]);
camlight;

left_surf_sm = smooth_cortex_fast(left_surf,.1,3000);
right_surf_sm = smooth_cortex_fast(right_surf,.1,3000);


figure;
patch('faces',left_surf_sm.faces,'Vertices',left_surf_sm.vertices,'facevertexcdata',a.datal,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(-90,0);clim([0,3]);
camlight;

figure;
patch('faces',right_surf_sm.faces,'Vertices',right_surf_sm.vertices,'facevertexcdata',a.datar,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(90,0);clim([0,3]);
camlight;

[subdir,q] = fileparts(mri_basefile);



% Activation on atlas (common template)
left_atlas_surf = readdfs([subdir, '/atlas.left.mid.cortex.svreg.dfs']);
right_atlas_surf = readdfs([subdir, '/atlas.right.mid.cortex.svreg.dfs']);


figure;
patch('faces',left_atlas_surf.faces,'Vertices',left_atlas_surf.vertices,'facevertexcdata',a.datal_atlas,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(-90,0);clim([0,3]);
camlight;

figure;
patch('faces',right_atlas_surf.faces,'Vertices',right_atlas_surf.vertices,'facevertexcdata',a.datar_atlas,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(90,0);clim([0,3]);
camlight;


left_atlas_surf_sm = smooth_cortex_fast(left_atlas_surf,.1,3000);
right_atlas_surf_sm = smooth_cortex_fast(right_atlas_surf,.1,3000);

figure;
patch('faces',left_atlas_surf_sm.faces,'Vertices',left_atlas_surf_sm.vertices,'facevertexcdata',a.datal_atlas,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(-90,0);clim([0,3]);
camlight;

figure;
patch('faces',right_atlas_surf_sm.faces,'Vertices',right_atlas_surf_sm.vertices,'facevertexcdata',a.datar_atlas,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(90,0);clim([0,3]);
camlight;


datal_atlas_sm = smooth_surf_function(left_atlas_surf, a.datal_atlas);
figure;
patch('faces',left_atlas_surf_sm.faces,'Vertices',left_atlas_surf_sm.vertices,'facevertexcdata',datal_atlas_sm,'edgecolor','none','facecolor','interp');
material dull;axis off;axis equal;axis tight;colormap('hot');view(-90,0);clim([0,3]);
camlight;

