# HTDPA Analyzer

**[Open the HTDA demonstration video](media/htda-demo_cmprsd.mp4)**

<video autoplay muted loop playsinline controls width="100%" poster="media/logo.png">
  <source src="./media/htda-demo_cmprsd.mp4" type="video/mp4">
  <a href="./media/htda-demo_cmprsd.mp4">Open the HTDA demonstration video</a>
</video>

<p align="center">
  <img src="media/logo.png" alt="Self-Driving Laboratory, Acceleration Consortium, University of Toronto" width="420">
</p>

**Self-Driving Laboratory (SDL5) Formulation**  
**Acceleration Consortium - University of Toronto**

This repository contains the data, software, device designs, laboratory automation, and computational analyses supporting a scientific manuscript submitted to *Digital Discovery*. It is organized as a data-availability and reproducibility package for the high-throughput diffusion/permeation apparatus (HTDPA) and its self-driving laboratory workflow.

> **Video behavior:** The demonstration is placed first and requests autoplay, mute, looping, and inline playback. Browsers and services such as GitHub may block autoplay or sanitize the HTML video element; the linked local MP4 remains available as a fallback. Opening a repository does not universally force a README to open or a video to play, so playback still depends on the host and browser.

## Repository map

| Path | Contents |
| --- | --- |
| [`Example_data/`](Example_data/) | Raw example release datasets in Excel format for formulation and probe conditions. |
| [`Programming/OT_Protocol/`](Programming/OT_Protocol/) | Opentrons protocol used for automated sampling, media replacement, and lid handling. |
| [`Programming/OT_LabWare/`](Programming/OT_LabWare/) | Custom labware definitions used by the OT-2 workflow. |
| [`CAD_files/`](CAD_files/) | STL geometry for the HTDPA cell, reciprocating pump, holders, and related components. |
| [`UI/`](UI/) | Streamlit dashboard for release-profile QC, kinetic fitting, and manuscript summary figures. |
| [`HTDA_Cell_CFD_Sim/`](HTDA_Cell_CFD_Sim/) | Transient Warp simulation, vorticity sweep, post-processing, notebook, and Ansys Workbench/Fluent files. |
| [`media/`](media/) | Repository presentation assets: the HTDA demonstration video and SDL5/Acceleration Consortium logo. |

## Reproducibility workflow

The repository is intended to preserve the path from physical design and automated experimentation to analysis and simulation:

1. **Device design:** Inspect the STL files in `CAD_files/` and the CFD/Workbench project in `HTDA_Cell_CFD_Sim/Workbench/`.
2. **Automated experiment:** Review the protocol and labware definitions in `Programming/` before loading a protocol onto an Opentrons system.
3. **Raw data:** Use the workbooks in `Example_data/` as the input data for release analysis. Preserve the original files when creating derived outputs.
4. **Analysis:** Run the dashboard in `UI/` to inspect release profiles, assay QC, kinetic models, apparent permeability, and the Korsmeyer-Peppas exponent.
5. **Flow simulation:** Follow [`HTDA_Cell_CFD_Sim/README.md`](HTDA_Cell_CFD_Sim/README.md) to reproduce the transient and steady-state design-supporting calculations.

## Release-analysis dashboard

The dashboard applies the publication-matched analysis rules documented in [`UI/README_HTDPA_dashboard.md`](UI/README_HTDPA_dashboard.md), including:

- strict-majority LOQ validation for quantitative time points;
- exclusion of negative cumulative-release values;
- common analysis windows for Zero-order, Higuchi, and First-order fits;
- optional start/end control for Korsmeyer-Peppas and Peppas-Sahlin fits;
- replicate-based CV reporting;
- six-panel summary output for release, variability, linearity, model fit, apparent `Kp`, and apparent `n`.

Install the dashboard dependencies and launch it from the dashboard directory:

```powershell
cd UI
python -m pip install -r requirements_htdpa_dashboard.txt
streamlit run htdpa_release_dashboard.py
```

The dashboard accepts the example Excel workbooks and can export analysis summaries as CSV and figures as high-resolution PNG files.

## CFD and flow-design analysis

The CFD subproject includes a GPU-capable NVIDIA Warp implementation of a transient 2D reciprocating-flow model, a flow-rate sweep, post-processing scripts, a Jupyter notebook, cached numerical outputs, and an Ansys Workbench/Fluent 3D design study. The precomputed outputs are included so the results can be inspected without rerunning the simulations.

```powershell
cd HTDA_Cell_CFD_Sim
conda env create -f environment.yml
conda activate htda-cfd
python run_all.py
```

See [`HTDA_Cell_CFD_Sim/README.md`](HTDA_Cell_CFD_Sim/README.md) for model parameters, individual stages, hardware requirements, and limitations.

## Data availability and provenance

- Raw example workbooks are retained in `Example_data/`.
- Analysis and simulation scripts are retained alongside their cached outputs and notebooks.
- CAD and mesh/project files are retained in `CAD_files/` and `HTDA_Cell_CFD_Sim/Workbench/`.
- Protocol and labware definitions are retained in `Programming/`.
- The dashboard changelog is maintained in [`UI/CHANGELOG.md`](UI/CHANGELOG.md).
- Generated figures and exported summaries should be treated as derived data and should retain a clear link to their input workbook and analysis settings.

This repository is a research artifact for the associated manuscript. Experimental interpretation should use the manuscript and the documented analysis assumptions together with the files here. The CFD model is design-supporting and includes the limitations described in its subproject README; it is not a substitute for experimental validation.

## Citation

Please cite the associated *Digital Discovery* manuscript and identify the repository revision used for analysis. A formal citation record can be added here when the manuscript and repository release receive their final bibliographic information.

## Contact

**Self-Driving Laboratory (SDL5) Formulation**  
Acceleration Consortium  
University of Toronto