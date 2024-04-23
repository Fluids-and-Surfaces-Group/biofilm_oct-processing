# Align the channel with the horizontal plane. Angles have to 
# be determined in advance.


from ij import IJ
#
from datetime import datetime
#
import os
#
# functions
#

def log(info = "", *additional):
	IJ.log(datetime.now().strftime('%d.%m.%Y, %H:%M:%S'))
	IJ.log(str(info))
	for i in additional:
		IJ.log(i)
#
def save_log(path, oct_file):
	filename, file_extension = os.path.splitext(oct_file)
	IJ.selectWindow("Log")
	_save = os.path.join(path, filename)
	IJ.saveAs("text", _save)
	IJ.run("Close")

#
# main program
#
path = IJ.getDirectory("Choose directory containing files...")
#
log("path: " + path, "\n")
# get file list
nFiles = 0

case = 2

if case==2:
# Iteration2
# Numbers determined manually
    angles_all = [[-20.8, -0.9], [-20.4, 1.1], [-20.4, 0], [-18.7, 0.7],
                    [-20, 0.2], [-20.5, -0.4], [-19.7, 0.5], [-19, -0.4], 
                    [-19.7, 0], [-20.3, 0.9],[-19.3, 0],[-18.2, 0.2]]

if case==3:
# Iteration3
    angles_all = [[-20.7, 0.9], [-20.7, 0.7], [-19.8, 1], [-18.4, 0], 
                    [-20.4, 0], [-20.7, -0.1], [-18.8, 1.2], [-17.6, 1.3],
                    [-20.4, -0.1], [-20.3, 1.1], [-19.5, 0.7], [-18.2, 0.9]]


# For all *.tif files in path
for _file in os.listdir(path):
    if _file.endswith(".tif"):
        # Get the file number
        fnumber = int(_file[8:11]) # file number
        
        
        if fnumber>=558 and fnumber <=575:
            # Get rotation angles for this channel
            angles = angles_all[nFiles%12]
            # Open file
            raw_file = _file
            log("Opening file: " + raw_file, "\n")
            IJ.open(os.path.join(path,raw_file))
            log("Opened file...")
            
            # Between file number 491 and 575, the scan orientation was different. To correct this, rotate scan by 90 degrees.
            if fnumber>=491 and fnumber <=575:
                # Correct orientation of faulty scans
                IJ.run("Reslice [/]...", "output=0.014 start=Left avoid");
                IJ.run("Rotate 90 Degrees Right");
            # Rotate first angle
            IJ.run("Rotate... ", "angle=" + str(angles[0]) + " grid=1 interpolation=Bilinear enlarge stack")
            # Between file number 491 and 575, the lateral resolution was different. Therefore, the pixel spacing is adjusted in this case.
            if fnumber>=491 and fnumber <=575:
                IJ.run("Reslice [/]...", "output=0.014 start=Left avoid")
            else:
                IJ.run("Reslice [/]...", "output=0.012 start=Left avoid")
            # Rotate second angle
            IJ.run("Rotate... ", "angle=" + str(angles[1]) + " grid=21 interpolation=Bilinear enlarge stack")
            if fnumber>=491 and fnumber <=575:
                IJ.run("Reslice [/]...", "output=0.014 start=Left avoid")
            else:
                IJ.run("Reslice [/]...", "output=0.012 start=Left avoid")
            # Save final image
            log("Saving as " + os.path.join(path, "rotated",_file))
            IJ.saveAs("Tiff", os.path.join(path, "rotated",_file))
            # Close all opened windows
            IJ.run("Close")
            IJ.run("Close")
            IJ.run("Close")
            if fnumber>=491 and fnumber <=575:
                # Close additional reslice
                IJ.run("Close")

        # raise counter
        nFiles += 1
#
IJ.showMessage("FINISHED", str(nFiles) + " OCT file(s) have successfully been processed!") 
