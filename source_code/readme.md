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

### 1. Set Parameter
* Click **Select Data Directory** → choose your main data folder (*Your_Data_Directory*)
      - Note: the main folder should contain subfolders with the name in the set: {Control, NTG, OAG} 
* In the **Disease** dropdown, select *CONTROL*, *NTG*, or *OAG*  (Currently only these three types are supported)
* In the **Subject ID** box, enter the subject folder name (e.g., *CV001*)
* In the **Timepoint(s)** box, enter subfolder names (e.g., *OD_MAC_08121030AM*)
* Review and adjust default parameters as needed  

Click **`Set Parameters`** (yellow button)
* Check Runtime Log: it should list the detected files

### 2. Multi-volume Registration
#### 2.1. Multi-volume registration
Click **`Multi-volume Registration`** (green button)
* In the pop-up Choose Template window
	- Select template files for each time point
	- Ensure en-face and volume (*_Cube*) templates are matched and correctly labeled  
	- A pop-up window will display registration progress
#### 2.2. 3D CoV calculation
* In the pop-up ROI input dialogue

	- Open *Subject_ID → Results3D → Timepoint/Average → \*__Angiography_3x3_FlowCube_z_all.tif* in ImageJ
	- Use rectangle tool to select ROI
	- Right-click rectangle → ROI Properties → check “List coordinates” → click OK
	- Copy the middle two Y-values into the ROI dialog in RetCoV GUI
	- A new window will pop up showing the progress of B-scan-based CoV computation  

**Outside GUI**
---------------
### 3. Retinal layer segmentation
* Layer boundaries from the innermost to the outermost need to be labelled with a number from 1-5:
   
     * ILM        = 1
     * NFL - GCL  = 2
     * IPL - INL  = 3
     * OPL - ONL  = 4
     
* Save your segmentation result as "Boundary.nii.gz"

### 4. Retinal Vessel Segmentation
* Vessel segmentation results need to be saved under folder:
  
	*SubjectID → Results3D → Timepoint → nifti*
* Please use the *\*_enfaceOCTA_fromZeiss.nii* as a reference when deriving the segmentation

**Go Back to GUI**
---------------
### 5. Generate Projection Images

Click **`Generate projection images`**

#### 5.1. Pop-up window: offset value
* Default can start with any number
	- Go back AFTER segmentation is done to check if the values are good, rerun this step and adjust the offset value if necessary)
 	- **How To Check**
		- Open *Results3D → timepoint → CoV → NFL_ONL → Binarized_avgOCTA_proj_offset_X.tif*
		- Drag the *.tif* file into ImageJ to view, compare with *\*_enfaceOCTA_fromVolume.tif* saved under the same folder
#### 5.2. Pop-up window: segmentation
* Select *Boundary.nii.gz* generated in **Step 3**


### 6. Threshold & visualization

Click **`Threshold & visualization`** to see histogram
- CoV values from instrument-processed data - choose tail point
- CoV values from volume-extracted data (choose same point as previous)
- CoV values from B-scan basic CoV computation (choose same point as previous)
- CoV values from B-scan basis CoV computation (deep vessels) - choose the intersection point (green circle)




Support
-------
For questions, bug reports, or feedback, please contact:
Yudan Chen  
University of British Columbia  
Email: cydan199@student.ubc.ca


