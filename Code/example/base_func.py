"""
==========================================================
This module contains utility functions for the influenza ABM simulation
built on the Covasim framework.

Author: Liping Peng
Institution: The University of Hong Kong
Date: 9 September 2025
"""


import pandas as pd
import os


# Set working directory - same in run_Season.py
path = "D:/ABM influenza/"
os.chdir(path)



def func_create_folder(path, folder_name):
    """
    Create a directory for storing simulation results if it doesn't exist.

    Args:
        path (str): Base path where the folder should be created
        folder_name (str): Name of the subfolder to create

    Returns:
        None: Creates directory structure on filesystem
    """
    folder_name = path + "/" + folder_name
    folder = os.path.exists(folder_name)
    if not folder:
        os.makedirs(folder_name)



def func_vac_distr(Season = None,  age_group = None):
    """
    Retrieve baseline vaccination coverage for specific season and age group.

        Args:
            Season (int): Influenza season number (1-6)
            age_group (str): Age group identifier (e.g., '0-11', '12-17', '18-49', etc.)

        Returns:
            float: Baseline vaccination coverage proportion for the specified group
    """
    data_vac = pd.read_csv("LP_input data/vac_coverage.csv")
    df = data_vac[data_vac['season'] == Season]
    vprop_base = df[age_group].values

    return vprop_base[0]

