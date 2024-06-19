% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.


function [nii] = load_nii_BIG_Lab(fname)
nii.img = niftiread(fname);

v_hdr = niftiinfo(fname);
nii.hdr = v_hdr;
nii.hdr.dime.pixdim = nii.hdr.raw.pixdim;
nii.untouch=1;
% nii.hdr.hist.magic = nii.hdr.magic;
% nii.hdr.dime.datatype = nii.hdr.datatype;
% nii.hdr.hk.sizeof_hdr = nii.hdr.sizeof_hdr;
% nii.hdr.hk.data_type = nii.hdr.datatype;
nii.hdr.hist.sform_code = nii.hdr.raw.sform_code;
nii.hdr.hist.qform_code =  nii.hdr.raw.qform_code;
nii.hdr.hist.srow_x= nii.hdr.raw.srow_x;
nii.hdr.hist.srow_y= nii.hdr.raw.srow_y;
nii.hdr.hist.srow_z= nii.hdr.raw.srow_z;
nii.hdr.hist.qoffset_x = nii.hdr.raw.qoffset_x;
nii.hdr.hist.qoffset_y = nii.hdr.raw.qoffset_y;
nii.hdr.hist.qoffset_z = nii.hdr.raw.qoffset_z;
nii.hdr.hist.quatern_b = nii.hdr.raw.quatern_b;
nii.hdr.hist.quatern_c = nii.hdr.raw.quatern_c;
nii.hdr.hist.quatern_d = nii.hdr.raw.quatern_d;

