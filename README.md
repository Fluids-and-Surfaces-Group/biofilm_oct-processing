# biofilm_oct-processing
These files are used to process OCT-scans, calculate statistical 
properties of the contained biofilm, and to plot the results.

# File processing
These scripts assume that the .tiff image stacks have been extracted from the .oct archives.
The following scripts should be run in sequence.

RotateOCTScans.py uses ImageJ to rotate the stacks along two axes, aligning the substratum with the horizontal plane. The angles need to be determined manually for each channel.

process_OCT_files.py detects the substratum, removes remaining warping due to optical effects, and binarises the image.

rescale_pixel_resolution.py adjusts a set of scans that were taken at a different lateral resolution. These scans are interpolated onto the resolution of the other measurements. 

register_OCT_scans.py is used to align scans of the same channel onto the same coordinate system, retaining only the segment that is present in all scans.

# Biofilm statistics
calculateStatisticsOCT.py calculates many biofilm parameters that are defined in utilities/biofilm_stats.py. These measures are then saved to a dataFrame.

The scripts beginning with plot_* create figures that are used in the paper using the files that were created in the prior steps.



# Copyright notice
Copyright 2024 Cornelius Wittig

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
