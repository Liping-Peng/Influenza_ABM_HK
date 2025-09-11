# Covasim Flu Modifications

This repository provides **modifications to a subset of Covasim (v3.1.4)** to enable **seasonal influenza modeling** with age-stratified transmission, improved immunity mechanisms, and calibration against age-specific incidence and school absenteeism data.

> **Note:** Only a subset of files are modified. Users should first install the original Covasim package and then replace the corresponding files with these modified versions.

---

## Overview

- Purpose: Adapt Covasim for influenza modeling, including age-specific transmission and calibration.
- Scope: Modifications affect only six files (see below). Other files remain under the original Covasim license.
- Transparency: All modifications are documented with the `LP_mod` prefix to ensure reproducibility.

---

## Modified Files

- `parameters.py`
- `analysis.py`
- `immunity.py`
- `people.py`
- `population.py`
- `sim.py`

All modifications are documented with the `LP_mod` prefix. See [MODIFICATIONS.md](MODIFICATIONS.md) for detailed descriptions of changes.

---

## Installation & Usage

1. **Install original Covasim (v3.1.4)**

```bash
git clone https://github.com/InstituteforDiseaseModeling/covasim.git
cd covasim
pip install -e .
```

2. **Apply the modifications**

Download the modified files from this repository:

* `parameters.py`
* `analysis.py`
* `immunity.py`
* `people.py`
* `population.py`
* `sim.py`

Replace the corresponding files in the original Covasim installation directory with these modified versions.

---

## License

* Original Covasim code remains under **IDM license** (© 2020–2022 Institute for Disease Modeling).
* Modifications in this repository are authored by **Liping Peng** and are shared under the ​**MIT License**​.

---

## Citation

If you use this repository in research, please cite both the original Covasim and this repository:

1. Original Covasim:
   Covasim: an agent-based model of COVID-19 dynamics and interventions. Kerr CC, Stuart RM, Mistry D, Abeysuriya RG, Rosenfeld R, Hart G, Núñez RC, Cohen JA, Selvaraj P, Hagedorn B, George L, Jastrzębski M, Izzo A, Fowler G, Palmer A, Delport D, Scott N, Kelly S, Bennette C, Wagner B, Chang S, Oron AP, Wenger E, Panovska-Griffiths J, Famulare M, Klein DJ (2021). *PLOS Computational Biology*​**17** (7): e1009149. doi: [https://doi.org/10.1371/journal.pcbi.1009149](https://doi.org/10.1371/journal.pcbi.1009149).
2. This repository:
   Peng L, Guo Y, Tsang NNY, Huang X, Hau C, Cowling BJ, Ip DKM, Tsang TK. Multi-source agent-based modeling to optimize influenza mitigation strategies in Hong Kong.

---

## Acknowledgements

This work builds upon ​**Covasim**​, developed by the Institute for Disease Modeling (IDM). All credit for the original model goes to the Covasim development team.

