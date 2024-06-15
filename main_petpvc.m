% Read the tissue fraction NIfTI files

tissueFracFile = '/home/ajoshi/Projects/PETBrainSuite/test_data/mri.pvc.frac.nii.gz'


tissueFracData = niftiread(tissueFracFile);

grayMatterData = (3-tissueFracData).*(tissueFracData>2) +(tissueFracData-1).*(tissueFracData<2).*(tissueFracData>=1);
whiteMatterData = (tissueFracData-2).*(tissueFracData>2);
CSFData = (2-tissueFracData).*(tissueFracData<=2).*(tissueFracData>1) + (tissueFracData).*(tissueFracData<=1);
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
niftiwrite(outputData, 'output_4D_tissue_fraction.nii.gz', niftiInfo);
