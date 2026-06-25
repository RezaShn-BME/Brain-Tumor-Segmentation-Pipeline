# Automated Brain Tumor Segmentation & Quantitative Reporting Pipeline

An end-to-end Biomedical Image Processing pipeline developed in MATLAB for automated head detection, brain masking, multi-focal tumor segmentation, and clinical-grade PDF report generation. This framework is specifically optimized to mitigate false positives in edge-case MRI scenarios (e.g., ring-enhancements, tilted slices, and post-operative healthy tissues).
---
## Key Engineering Features

**Adaptive Brain Masking:** Implements multi-level thresholding and morphological operations (`imerode`, `imfill`) to isolate the brain zone, ensuring robustness against tilted imaging planes.[cite: 2]
**Pathophysiological Scoring System:** Employs a dual-metric scoring system based on cluster Area x Solidity to accurately capture multi-focal tumor regions while filtering out high-frequency noise.[cite: 2]
**Advanced Intensity Checking:** Incorporates background brain statistical variance (Mean & Standard Deviation) to protect healthy/post-op slices from false-positive segmentation.[cite: 2]
**Quantitative Research Validation:** Integrates automatic validation using **Dice Similarity Coefficient** and **Jaccard Index (IoU)** against Ground Truth masks (Optional).[cite: 2]
**Automated Clinical Reporting:** Dynamically generates an official PDF analytics report utilizing the MATLAB Report Generator API (`mlreportgen`).[cite: 2]
---
## Repository Structure
```text
├── src/
│   └── detectAndReportTumor.m    # Core MATLAB pipeline function
├── data/
│   └── sample_mri.png            # Example input brain MRI slice
└── README.md                     # Project documentation
```
---
## Pipeline Architecture & Algorithmic Steps

* **Preprocessing:** Grayscale conversion and Adaptive Median Filtering (`medfilt2`) for salt-and-pepper noise reduction.

* **Tissue Masking:** Head segmentation via global Otsu's thresholding, followed by dynamic bounding-box cropping to define the brain-only region of interest (ROI).

* **Adaptive Segmentation:** Multi-level thresholding (`multithresh`) coupled with holistic hole-filling to capture necrotic cores in ring-enhancing lesions.

* **Statistical Validation:** Regional component evaluation compared against normal brain tissue parameters (Mean + 1.2 x STD).

* **Report Generation:** Structural exporting of volumetric percentages, bounding-box parameters, and spatial overlays into a formatted PDF.

---

## Usage & Execution

To execute the pipeline with validation metrics, run the following command in the MATLAB Command Window:

```matlab
% Define paths
imagePath = 'data/sample_mri.png';
authorName = 'Reza Shojaei Nasab';

% Run Pipeline
detectAndReportTumor(imagePath, authorName);


## Contact & Research Collaborations

* **Developer:** Reza Shojaei Nasab
* **Field:** Biomedical Engineering (Bioelectric)
---
## Results Visualization
![Tumor Detection Result](output/Tumor_Result.png)
