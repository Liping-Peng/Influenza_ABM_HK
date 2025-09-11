# **Influenza_ABM_HK**

Influenza_ABM_HK is an agent-based modeling framework adapted from Covasim (v3.1.4) to simulate **seasonal influenza dynamics** across six Hong Kong seasons. It calibrated using ​**age-stratified incidence and school absenteeism data**​. The framework incorporates ​age-specific transmission, immunity mechanisms, and subgroup-targeted interventions​, enabling realistic epidemic simulations and evaluation of ​**age-targeted vaccination programs, stay-at-home policies, and  mask-wearing**​.


---

## **Repository Structure**

```text
respiratory/
├── code/
│   ├── Covasim_Flu_Modification/  # Modified Covasim framework for influenza
│   │   ├── README.md              # Details of modifications
│   │   ├── MODIFICATIONS.md       # Detailed change log
│   │   ├── parameters.py
│   │   ├── analysis.py
│   │   ├── immunity.py
│   │   ├── people.py
│   │   ├── population.py
│   │   └── sim.py
│   └── example/                   # Example scripts for running simulations
│       ├── README.md              # Instructions for example usage
│       ├── run_Season.py          # Main simulation script
│       ├── base_func.py           # Utility functions
│       ├── base_data.py           # Seasonal datasets
│       └── base_subgroup.py       # Population subgroup definitions for interventions
└── data/                          # Required input datasets
```
---

## **Features**

- Influenza-specific modifications to disease progression and prognoses
- Age-stratified transmission dynamics
- Calibration using Hong Kong age-specific incidence and school absenteeism data
- Flexible vaccine types for both match and mismatch seasons
- Subgroup-targeted interventions to evaluate population-specific strategies

---

## **Getting Started**

### 1. Install Covasim

Install the original Covasim package (v3.1.4):

```bash
git clone https://github.com/InstituteforDiseaseModeling/covasim.git
cd covasim
pip install -e .
```

### 2. Apply Modifications

Replace the following core files in the original Covasim installation with the modified versions:

- `parameters.py`
- `analysis.py`
- `immunity.py`
- `people.py`
- `population.py`
- `sim.py`

### 3. Run Example Simulations

Navigate to the example folder:

```
cd code/example
python run_Season.py
```

Modify parameters, scenarios, or interventions in run_Season.py as needed

---

## **License**

### Original Covasim:

- Covasim 3.1.4 (2022-10-22) — © 2020–2022 Institute for Disease Modeling (IDM).

- All unmodified files remain under the original Covasim license.

### Modifications in this repository:

- Only a subset of files are modified: `parameters.py`, `analysis.py`, `immunity.py`, `people.py`, `population.py`, and `sim.py`.
- Modifications are documented with the `LP_mod` prefix for transparency.

- These modifications are authored by Liping Peng and may be cited in academic work.

---

## **Citation**


If you use this repository in research, please cite both the original Covasim and this repository:

#### 1. Original Covasim:

Covasim: an agent-based model of COVID-19 dynamics and interventions. Kerr CC, Stuart RM, Mistry D, Abeysuriya RG, Rosenfeld R, Hart G, Núñez RC, Cohen JA, Selvaraj P, Hagedorn B, George L, Jastrzębski M, Izzo A, Fowler G, Palmer A, Delport D, Scott N, Kelly S, Bennette C, Wagner B, Chang S, Oron AP, Wenger E, Panovska-Griffiths J, Famulare M, Klein DJ (2021). *PLOS Computational Biology*​**17** (7): e1009149. doi: [https://doi.org/10.1371/journal.pcbi.1009149](https://doi.org/10.1371/journal.pcbi.1009149).

#### 2. This repository:

Peng L, Guo Y, Tsang NNY, Huang X, Hau C, Cowling BJ, Ip DKM, Tsang TK. Multi-source agent-based modeling to optimize influenza mitigation strategies in Hong Kong.

---

## **Acknowledgements**

This work builds upon ​**Covasim**​, developed by the Institute for Disease Modeling (IDM). All credit for the original model goes to the Covasim development team.



