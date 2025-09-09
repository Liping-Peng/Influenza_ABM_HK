"""
==========================================================================
This module contains subtargeting functions for prioritizing different interventions
across population subgroups. These functions enable age-specific targeting of
vaccination and stay-at-home interventions.

Author: Liping Peng
Institution: The University of Hong Kong
Date: 9 September 2025
"""

import numpy as np
import covasim as cv
import base_func as bfc



def Sub_vac_cover(sim):
    """
    Subtargeting function for vaccination coverage by age group and strategy.

    Args:
        sim (cv.Sim): The current simulation object

    Returns:
        dict: Dictionary with indices and vaccination probability values for each agent

    Note:
        Implements different vaccination strategies (vg0-vg7) with age-specific coverage
        based on historical data or intervention parameters.
    """

    # Get vaccination strategy parameters from simulation
    Season = sim.people.season
    param_vg = sim.pars['LP_param_vg']
    param_vc = sim.pars['LP_param_vc']

    # Define age groups for vaccination targeting
    age_0_5 = cv.true(sim.people.age < 6)
    age_6_11 = cv.true((sim.people.age >= 6) * (sim.people.age < 12))
    age_12_17 = cv.true((sim.people.age >= 12) * (sim.people.age < 18))
    age_18_39 = cv.true((sim.people.age >= 18) * (sim.people.age < 40))
    age_40_64 = cv.true((sim.people.age >= 40) * (sim.people.age < 65))
    age_65 = cv.true(sim.people.age >= 65)
    inds = sim.people.uid  # Everyone in the population -- equivalent to np.arange(len(sim.people))
    vals = np.ones(len(sim.people))  # Create the array

    # Set baseline vaccination coverage from historical data
    vals[age_0_5] = bfc.func_vac_distr(Season = Season, age_group = '0_5')
    vals[age_6_11] = bfc.func_vac_distr(Season = Season,  age_group = '6_11')
    vals[age_12_17] = bfc.func_vac_distr(Season = Season,  age_group = '12_17')
    vals[age_18_39] = bfc.func_vac_distr(Season = Season,  age_group = '18_39')
    vals[age_40_64] = bfc.func_vac_distr(Season = Season,  age_group = '40_64')
    vals[age_65] = bfc.func_vac_distr(Season = Season,  age_group = '>=65')

    # Apply vaccination strategy interventions
    if param_vg == 'vg0':
        # Baseline strategy: no additional coverage
        pass
    elif param_vg == 'vg1':
        # Strategy 1: Expand to children under 12
        vals[age_0_5] = bfc.func_vac_distr(Season = Season,  age_group = '0_5') + param_vc
        vals[age_6_11] = bfc.func_vac_distr(Season = Season,  age_group = '6_11') + param_vc
    elif param_vg == 'vg2':
        # Strategy 2: Expand to children under 18
        vals[age_0_5] = bfc.func_vac_distr(Season = Season,  age_group = '0_5') + param_vc
        vals[age_6_11] = bfc.func_vac_distr(Season = Season,  age_group = '6_11') + param_vc
        vals[age_12_17] = bfc.func_vac_distr(Season = Season,  age_group = '12_17') + param_vc
    elif param_vg == 'vg3':
        # Strategy 3: Expand to seniors 65+
        vals[age_65] = bfc.func_vac_distr(Season = Season,  age_group = '>=65') + param_vc
    elif param_vg == 'vg4':
        # Strategy 4: Expand to children and seniors
        vals[age_0_5] = bfc.func_vac_distr(Season = Season,  age_group = '0_5') + param_vc
        vals[age_6_11] = bfc.func_vac_distr(Season = Season,  age_group = '6_11') + param_vc
        vals[age_12_17] = bfc.func_vac_distr(Season = Season,  age_group = '12_17') + param_vc
        vals[age_65] = bfc.func_vac_distr(Season = Season,  age_group = '>=65') + param_vc
    elif param_vg == 'vg5':
        # Strategy 5: Universal vaccination
        vals[age_0_5] = bfc.func_vac_distr(Season = Season,  age_group = '0_5') + param_vc
        vals[age_6_11] = bfc.func_vac_distr(Season = Season,  age_group = '6_11') + param_vc
        vals[age_12_17] = bfc.func_vac_distr(Season = Season,  age_group = '12_17') + param_vc
        vals[age_18_39] = bfc.func_vac_distr(Season = Season,  age_group = '18_39') + param_vc
        vals[age_40_64] = bfc.func_vac_distr(Season = Season,  age_group = '40_64') + param_vc
        vals[age_65] = bfc.func_vac_distr(Season = Season,  age_group = '>=65') + param_vc
    elif param_vg == 'vg6':
        # Strategy 6: Use 2024 vaccination coverage rates
        vals[age_0_5] = bfc.func_vac_distr(Season = 2024,  age_group = '0_5')
        vals[age_6_11] = bfc.func_vac_distr(Season = 2024,  age_group = '6_11')
        vals[age_12_17] = bfc.func_vac_distr(Season = 2024,  age_group = '12_17')
        vals[age_18_39] = bfc.func_vac_distr(Season = 2024,  age_group = '18_39')
        vals[age_40_64] = bfc.func_vac_distr(Season = 2024,  age_group = '40_64')
        vals[age_65] = bfc.func_vac_distr(Season = 2024,  age_group = '>=65')
    elif param_vg == 'vg7':
        # Strategy 7: Target school-age children (3-17 years)
        age_3_17 = cv.true((sim.people.age >= 3) * (sim.people.age < 18))
        # Option 1: Add to baseline coverage
        # vals[age_3_17] = vals[age_3_17] + param_vc
        # Option 2: Set absolute coverage level (0-100%)
        vals[age_3_17] = param_vc
    else:
        raise ValueError("Vaccination strategy is not defined")


    # Return indices and values for vaccination intervention
    output = dict(inds=inds, vals=vals)
    return output





def Sub_tip_base(sim):
    """
    Subtargeting function for baseline stay-at-home probabilities.

    Args:
        sim (cv.Sim): The current simulation object

    Returns:
        dict: Dictionary with indices and probability values for each agent

    Note:
        Implements age-specific stay-at-home probabilities based on health-seeking behavior
        research (Zhang et al., BMC Public Health 2020).
    """

    # Define symptomatic age groups for stay-at-home targeting
    age_0_15 = cv.true((sim.people.age <= 15) * (sim.people.symptomatic))
    age_16_54 = cv.true((sim.people.age >= 16) * (sim.people.age <= 54) * (sim.people.symptomatic))
    age_55 = cv.true((sim.people.age >= 55) * (sim.people.symptomatic))
    inds = sim.people.uid  # Everyone in the population -- equivalent to np.arange(len(sim.people))
    vals = np.zeros(len(sim.people))  # Create the array

    # Set age-specific daily stay-at-home probabilities (ref: Q. Zhang.BMC PH.2020.)
    vals[age_0_15] = 0.181 * sim.pars['LP_base_tip']   # Daily prob： 1-（1-0.181)^6 = 0.699，1-（1-0.065)^6 = 0.332，1-（1-0.095)^6 = 0.449
    vals[age_16_54] = 0.065 * sim.pars['LP_base_tip']
    vals[age_55] = 0.095 * sim.pars['LP_base_tip']

    output = dict(inds=inds, vals=vals)
    return output





def Sub_tip_param(sim):
    """
    Subtargeting function for enhanced stay-at-home probabilities.

    Args:
        sim (cv.Sim): The current simulation object

    Returns:
        dict: Dictionary with indices and probability values for each agent

    Note:
        Extends baseline probability with intervention-specific improvements.
    """

    # Get intervention parameter for improvement
    param_tip = sim.pars['LP_param_tip']

    # Define symptomatic age groups for stay-at-home targeting
    age_0_15 = cv.true((sim.people.age <= 15) * (sim.people.symptomatic))
    age_16_54 = cv.true((sim.people.age >= 16) * (sim.people.age <= 54) * (sim.people.symptomatic))
    age_55 = cv.true((sim.people.age >= 55) * (sim.people.symptomatic))
    inds = sim.people.uid
    vals = np.zeros(len(sim.people))

    # Set enhanced age-specific daily probabilities
    # multiplied by (baseline + intervention improvement) to determine final daily stay-at-homes probability for symptomatic cases
    vals[age_0_15] = 0.181 * (sim.pars['LP_base_tip'] + param_tip)
    vals[age_16_54] = 0.065 * (sim.pars['LP_base_tip'] + param_tip)
    vals[age_55] = 0.095 * (sim.pars['LP_base_tip'] + param_tip)

    output = dict(inds=inds, vals=vals)
    return output