RetCOV User Manual (v1.0)
=========================

Overview
---------
RetCOV is a processing framework designed to perform 3D Coefficient of Variation (CoV) analysis for quantifying retinal capillary perfusion variability using OCT/OCTA imaging data.

System Requirements
--------------------------------
- Operating System: Windows 10 or later (64-bit)
- MATLAB Runtime: R2023b (automatically prompted during installation if not installed)

Required Toolboxes (for source users)
-------------------------------------------------------
Note: These toolboxes are only necessary if running the app from source code (.mlapp). They are **not** required for end users installing the standalone app.

- Simulink
- Deep Learning Toolbox
- Image Processing Toolbox
- Computer Vision Toolbox
- Natural-Order Filename Sort
- Curve Fitting Toolbox
- Econometrics Toolbox
- Sensor Fusion and Tracking Toolbox

Usage (for source users)
---------
 Double-click `RetCoV.mlapp` to open the GUI. Make sure the path of helper functions is added to the MATLAB path.

1. Set Parameter
   a. Click **Select Data Directory** → choose your main data folder (`Your_Data_Directory`)
      - Note: the main folder should contain subfolders with the name in the set: {Control, NTG, OAG} 
   b. In the ‘Disease’ dropdown, select `CONTROL`, `NTG`, or `OAG`  (Currently only these three types are supported)
   c. In the ‘Subject ID’ box, enter the subject folder name (e.g., `CV001`)  
   d. In the ‘Timepoint(s)’ box, enter subfolder names (e.g., `OD_MAC_08121030AM`)  
   e. Review and adjust default parameters as needed  
   Click **‘Set Parameters’** (yellow button)  
      - Check Runtime Log: it should list the detected files

2. Multi-volume Registration
   2.1. Multi-volume registration
	Click **‘Multi-volume Registration’** (green button)
	- In the pop-up menu, select template files for each time point  
   	   a. Ensure en-face and volume (_Cube) templates are matched and correctly labeled  
   	   b. A pop-up window will display registration progress
   2.2. 3D CoV calculation
	- ROI input dialog
	   a. Open `Subject_ID/Results3D/Timepoint/Average/*__Angiography_3x3_FlowCube_z_all.tif` in ImageJ  
	   b. Use rectangle tool to select ROI  
	   c. Right-click rectangle → ROI Properties → check “List coordinates” → click OK  
	   d. Copy the middle two Y-values into the ROI dialog in RetCoV GUI  
	   e. A new window will pop up showing the progress of B-scan-based CoV computation  


3. Retinal layer segmentation
   - layer boundary from the innermost to the outermost need to be label with number from 1-5:
   
     ILM        = 1
     NFL - GCL  = 2
     IPL - INL  = 3
     OPL - ONL  = 4
     
   - Save your segmentation result as "Boundary.nii.gz"

4. Retinal Vessel Segmentation
   a. vessel segmentation results need to be saved under folder:
      SubjectID → Results3D → Timepoint → nifty
   b. Please use the "*_enfaceOCTA_fromZeiss.nii" as a reference when deriving the segmentation

5. Generate Projection Images
   a. Click **Generate projection images**
   4.1. Pop-up window: offset value
        a. Default start with ‘10’ → go back AFTER segmentation is done to check if the values are good, rerun this step and adjust if necessary)
           ** To check **
              i. Results3D → timepoint → CoV → NFL_ONL → Binarized_avgOCTA_proj_offset_X.tif
              ii. Drag file into ImageJ to view, compare with “*_enfaceOCTA_fromVolume.tif”
   4.2. Pop-up window: segmentation
        a. select ‘Boundary.nii.gz’


6. Threshold & visualization
   a. Click ** Threshold & visualization ** to see histogram
CoV values from instrument-processed data - choose tail point
CoV values from volume-extracted data - choose same point as a)
CoV values from B-scan basic CoV computation -choose same point as a)
CoV values from B-scan basis CoV computation (deep vessels) - choose the intersection point (green circle)




Support
-------
For questions, bug reports, or feedback, please contact:
Yudan Chen  
University of British Columbia  
Email: cydan199@student.ubc.ca


