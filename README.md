# Hybrid ensemble Kalman--particle filtering for COVID-19 behavioral dynamics: Insights from South Korea, Japan, and Pakistan

This repository contains the code used to generate the figures and supporting analyses for the manuscript:

> **Hybrid ensemble Kalman--particle filtering for COVID-19 behavioral dynamics: Insights from South Korea, Japan, and Pakistan**

## Overview

Understanding the interplay between human behavioral responses and infectious disease transmission is essential for improving epidemic forecasting and evaluating intervention strategies. This study develops a hybrid ensemble Kalman--particle filtering (EnKPF) framework to estimate behavioral parameters and epidemic states from COVID-19 surveillance data.

The proposed approach combines the computational efficiency of the Ensemble Kalman Filter (EnKF) with the flexibility of Particle Filtering (PF) to capture nonlinear and non-Gaussian epidemic dynamics. The framework is applied to COVID-19 data from:

- South Korea
- Japan
- Pakistan

The study investigates:

- Time-varying behavioral transmission dynamics,
- Behavioral dissipation mechanisms,
- Intervention effectiveness across epidemic waves,
- Cross-country comparisons of behavioral responses, and
- Sensitivity of model outputs to behavioral parameters.

---

## Repository Structure

```text
.
├── manuscript_figures.executed.ipynb    # Notebook used to generate manuscript figures
├── data/                                # Input datasets and processed outputs
├── figures/                             # Generated manuscript figures
├── README.md
└── supplementary_files/                 # Additional supporting material (if applicable)
```

---

## Figures

The notebook

```text
manuscript_figures.executed.ipynb
```

reproduces the figures presented in the manuscript, including:

- Parameter evolution across epidemic waves,
- Behavioral transmission and dissipation estimates,
- Policy impact analyses,
- Intervention comparison heatmaps,
- Global sensitivity analyses using Partial Rank Correlation Coefficients (PRCC),
- Additional supplementary figures.

---

## Requirements

The analyses were performed using Python.

### Recommended Python Version

```text
Python 3.12+
```

### Main Dependencies

```text
numpy
scipy
pandas
matplotlib
seaborn
jupyter
```

Install the required packages using:

```bash
pip install numpy scipy pandas matplotlib seaborn jupyter
```

Alternatively:

```bash
pip install -r requirements.txt
```

if a `requirements.txt` file is provided.

---

## Usage

Clone the repository:

```bash
git clone https://github.com/your-username/your-repository.git
cd your-repository
```

Launch Jupyter Notebook:

```bash
jupyter notebook
```

Open and execute:

```text
manuscript_figures.executed.ipynb
```

The notebook will load the processed datasets, reproduce the analyses, and generate the manuscript figures.

---

## Data Availability

The repository includes the processed datasets required to reproduce the published figures. If restrictions apply to any raw data sources, users should obtain these datasets from the corresponding public health authorities or official repositories.

---

## Citation

If you use this code or build upon this work, please cite:

```text
Wasim Abbas et al.
Hybrid ensemble Kalman--particle filtering for COVID-19 behavioral dynamics:
Insights from South Korea, Japan, and Pakistan.
(Under review)
```




## Correspondence

**Dr. Wasim Abbas**  
Postdoctoral Researcher  
Kyungpook National University, Republic of Korea

