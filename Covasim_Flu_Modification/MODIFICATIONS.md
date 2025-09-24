## **Modified Files**

### **1. parameters.py**

#### 1.1New Parameters

|  **Parameter**                                                                                  |  **Description**                                                                                                      |
| -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
|  pars['LP\_season']                                                                              |  Influenza  season identifier (1–6).                                                                                  |
|  pars['LP\_flutype']                                                                             |  Dominant  circulating subtype in the simulation ('H1N1' or 'H3N2').                                                   |
|  pars['LP\_imprinting']                                                                          |  Boolean  flag for simulating immune imprinting effects.                                                               |
|  pars['LP\_preHAI']                                                                              |  Boolean  flag for pre-assigning HAI titers at initialization.                                                         |
|  pars['LP\_symp\_probs']                                                                         |  Custom  age-specific symptomatic probabilities (overrides defaults).                                                  |
|  pars['LP\_mismatch\_age']                                                                       |  Whether  to include age-specific incidence data for calibration.                                                      |
|  pars['LP\_mismatch\_absent']                                                                    |  Whether  to include school absenteeism data for calibration.                                                          |
|  pars['LP\_trans\_ORs']                                                                          |  Transmission  odds ratios by age group [0–12, 13–25, 26–45, 46–65, 65+]. Default = [2.0,  1.0, 1.0, 1.0, 1.0].    |
|  pars['LP\_sus\_ORs']                                                                            |  Susceptibility  odds ratios by age group [0–12, 13–25, 26–45, 46–65, 65+]. Default = [1.0,  1.0, 1.0, 1.0, 1.0].  |
|  pars['LP\_sus\_0\_25'],  pars['LP\_sus\_25\_45'], pars['LP\_sus\_45\_65'], pars['LP\_sus\_65']  |  Relative  susceptibility for coarse age bins.                                                                         |
|  pars['LP\_param\_vg']                                                                           |  Targeted  vaccination group (vg0: none, vg5: universal).                                                              |
|  pars['LP\_param\_vc']                                                                           |  Increase  in vaccine coverage (percentage points) above baseline.                                                     |
|  pars['LP\_base\_tip']                                                                           |  Baseline  daily probability of adherence to stay-at-home.                                                             |
|  pars['LP\_param\_tip']                                                                          |  Incremental  increase in stay-at-home adherence probability.                                                          |
|  pars['LP\_end\_tida']                                                                           |  Early  termination of stay-at-home (days before recovery).                                                            |

<pre><i><span><o:p> </o:p></span></i></pre>

#### 1.2 Modified Existing Parameters

- **beta_layer**: updated to `dict(h=2.0, s=0.5, w=0.5, c=0.2)`
- **contacts**: updated to `dict(h=2.0, s=18, w=18, c=20)`

#### 1.3 Complete Prognoses Overhaul

Replaced COVID-19 disease progression with influenza-specific parameters:

- Reduced age groups from 10 to 5 bins: `[0, 12, 25, 45, 65]`
- Set all severe outcomes to zero: `severe_probs = crit_probs = death_probs = 0.0`
- Focus on transmission dynamics only

#### 1.4 New Vaccine Types

Extended vaccine choices for influenza modeling:

- **lp_flu_30**: 30% efficacy vaccine for mismatch seasons
- **lp_flu_70**: 70% efficacy vaccine for match seasons

#### 1.5 New Helper Functions

Added three utility functions for advanced immune modeling:

- `get_imprt_data()`: Loads immune imprinting probability data
- `get_preHAI_data()`: Loads HAI titer distribution data
- `randomize_imprt()`: Stochastic assignment of imprinting status with ratio-based protection

##### 1.6 Data Dependencies

Requires external data files in `LP_input data/` directory:

- `imprinting_prob/sample_year_2019.csv`
- `hai_protection/HAI_distribution.csv`

```

```

### **2. analysis.py**

#### 2.1 Summary of Modifications

|  **Feature**                                      |  **Modification**                                                                                                                                                                                                                                                                           |
| ---------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|  Goodness-of-fit  (GOF) calculation                |  -  Explicitly pass the sim object into compute() and compute\_gofs(). -  Replaced normalized absolute error with Mean Squared Error (MSE) in compute\_gofs().                                                                                                                               |
|  Calibration  against age-specific incidence data  |  -  Added calculation of daily new infections for age groups (0–24, 25–44, 45–64,  ≥65) using infection\_log. -  Compute GOF metrics for each age group when sim.pars['LP\_mismatch\_age'] =  True. -  Added visualization of age-specific calibration results (calib.LP\_plot\_age()).  |
|  Calibration  against school absenteeism records   |  -  Implemented calculation of school absenteeism (ages 6–18) from custom  isolation logs (LP\_iso\_log). -  Compute GOF for absenteeism when sim.pars['LP\_mismatch\_absent'] = True. -  Added visualization of absenteeism calibration results (calib.LP\_plot\_absent()).                |
|  Configuration                                     |  -  Introduced control parameters: • sim.pars['LP\_mismatch\_age']  → enable age-specific fit • sim.pars['LP\_mismatch\_absent']  → enable absenteeism fit  • sim.pars['LP\_season'] = 1/2/3/4/5/6 →  specify influenza season for data handling                                       |

```

```

### **3. immunity.py**

#### 3.1 Summary of Modifications

|  **Feature**                |  **Modification**                                                                                                                                                                                                    |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|  Pre-existing  immunity      |  -  Added hemagglutination inhibition (HAI) titer-based immunity. -  Introduced discrete HAI titers (0–9 scale). -  Added preHAI\_bool parameter to enable/disable pre-assignment of HAI titers at  initialization.  |
|  Immune  imprinting          |  -  Implemented calc\_imprt() to simulate immune imprinting effects. -  Added imprt\_bool flag (currently inactive) for toggling imprinting.                                                                          |
|  Immunity  boosting          |  -  Replaced multiplicative boosting with an additive scheme: peak\_nab +  boost\_factor. -  Standardized to four-fold HAI titer rise upon re-exposure.                                                               |
|  Waning  pattern             |  -  Modified nab\_growth\_decay() function. -  Replaced exponential waning with linear decay at 14% per year.                                                                                                         |
|  Immunity  update functions  |  -  Overhauled immunity update mechanics: •  Modified update\_nab() •  Modified update\_peak\_nab()                                                                                                                 |

```

```

### **4. people.py**

#### 4.1 Summary of Modifications

Extended the infection logging system to capture detailed epidemiological metadata, enabling comprehensive reconstruction of transmission networks and supporting age-structured transmission dynamics analysis for advanced epidemiological research.

#### 4.2 New Infection Log Fields

|  **Field**         |  **Description**                     |
| --------------------- | --------------------------------------- |
|  age\_source        |  Age  of the **infector** individual   |
|  age\_target        |  Age  of the **infectee** individual   |
|  date\_symptomatic  |  Symptom  onset date of the infectee  |

```

```

### **5. population.py**

#### 5.1 Summary of Modifications

Extended the population initialization pipeline to propagate influenza-specific parameters throughout the model. Key function signatures were modified to accept influenza-related arguments, ensuring compatibility with age-structured influenza epidemiology and Hong Kong–specific synthetic population data.

#### 5.2   Functions Modified

##### (1) `make_people()`

- Added parameters: `LP_season`, `LP_sus_ORs`, `LP_trans_ORs`, `LP_symp_probs`.
- Default values for prognostic variables replaced with influenza-specific parameters:
  - `sim['prognoses']['sus_ORs'] = LP_sus_ORs`
  - `sim['prognoses']['symp_probs'] = LP_symp_probs`
  - `sim['prognoses']['trans_ORs'] = LP_trans_ORs`

##### (2) `make_synthpop()`

- Added `LP_season` parameter to enable season-specific synthetic population generation, tailored to Hong Kong demographics.

#### 5.3 SynthPops Population Generation Settings

|  **Parameter**                                    |  **Value**                                  |  **Description**                          |
| ---------------------------------------------------- | ---------------------------------------------- | -------------------------------------------- |
|  `with_school_types`              |  `True`                     |  Creates  explicit school types            |
|  `school_mixing_type`             |  `'age_clustered'`          |  Age-clustered  mixing within schools      |
|  `average_class_size`             |  `27`                       |  Based  on Hong Kong education statistics  |
|  `average_student_teacher_ratio`  |  `13`                       |  Composite  ratio derived from HK data     |
|  `country_location`               |  `'China%s'% LP_season`  |  Season-specific  age distributions        |
|  `state_location`                 |  `"Hongkong"`               |  Hong  Kong–specific configuration        |
|  `location`                       |  `"Hongkong"`               |  Regional  identifier for consistency      |

```

```

### **6. sim.py**

#### 6.1 Summary of Modifications

Enhanced seasonal influenza modeling with age-stratified transmission dynamics and updated immunity mechanisms, improving realism in population-level simulations.

#### 6.2   Modified Functions

##### (1) `init_people()`

- Incorporated age-stratified transmission and susceptibility parameters.
- Key changes:
  - Extracts user-defined parameters: `LP_sus_ORs`, `LP_trans_ORs`.
  - Passes age-specific parameters into `cvpop.make_people()`.
  - Stores influenza-specific immunity and intervention parameters (e.g., `imprt_bool`, `preHAI_bool`, `LP_season`) within `self.people` for downstream use.

##### (2) `step()`

- Introduced prevention of reinfection by providing complete immunity to previously exposed individuals.
- Key changes:
  - Identifies agents with `naive == False` (already exposed).
  - Sets their susceptibility immunity (`sus_imm`) to 1.0, corresponding to 100% protection.
    

