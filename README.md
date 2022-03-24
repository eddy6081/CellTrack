# CellTrack
Matlab and Python implementations of cell tracking from fluorescence microscopy images

## Matlab Implementation
Matlab implementation tracks cells using a simple minimum distance scheme. Scripts provided are meant to combine partial segments of cells in binary stacks, and then track and perform measures on full 2D projections of cells. Please see videos for usage!

**General workflow**:

(1) Use Overlap_v3.m to combine partial detections of cells in several manually segmented stacks. In anytime an intersection is detected between cells of consecutive image stacks, the user is requested to specify if they are one or two cells.
- Uses segmented binary images for z-stack images. Generates folders containing overlapped cells and non-overlapped cells. If there are several prepared folders you are interested to find overlapping cells in, follow this guide:
- folders: 1-3, 4-6, 7-9.
- Run Overlap_v3 on 1-3 and 4-6. This creates a folder 1-3_4-6_nonoverlap, which contains the cells in folder 1-3 that is not shared in 4-6. Additionally, folders 4-6_1-3_nonoverlap (same process) and 1-3_4_6_overlap (which contains cells that are SHARED between 1-3 and 4-6) are made.
- Next, run Overlap_v3 on 4-6_1-3_nonoverlap and 7-9. Again, it generates the same types of folders where 7-9_4-6_1-3_overlap has the overlapped cells between layers 4-6 and 7-9.

| Step 1 | Step 2 |
|---|---|
| ![](/Matlab/Tutorial/Overlap_1.mov) | ![](/Matlab/Tutorial/Overlap_2.mov) |


(2) Use PreviousCell_2D_v2.m to track and measure cell geometries.
- Requires maximum projection images and segmented cells. Recursively asks user to indicate which image stacks are to be analyzed. Asks user to select a folder in which to save cell-labeled images during tracking.
- If the user selects to screen cell detections, then program uses a guided user interface to allow user to select cells one wishes to follow or be analyzed in current and future frames.
- Returns saved .mat files in the same folder as the cell-labeled images.
- Once complete, returns a .csv file with cell labels and their geometric measurements.

![](/Matlab/Tutorial/PreviousCell_1.mov)

(3a) Use Detect_Error.m which uses cell centroid locations to make user aware of possible errors during tracking and request user to specify how to correct detections. Generates a .txt file that tracks any changes that were made, and creates a new .csv file with user changes.

![](/Matlab/Tutorial/Detect_Error.mov)

(3b) After reviewing the tracked cell-labeled images, if a mistake was made that was not detected by Detect_Error.m, or if a cell is mislabeled, use ChangeCellNum.m to fix the error.

(4) If step (3) was used, use Correct_CellNum_Ims.m to save edited cell-labeled images. Repeat steps 3 and 4 as necessary.
