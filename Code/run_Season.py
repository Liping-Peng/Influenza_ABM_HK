"""
Influenza ABM Simulation
==========================================================

This module implements an agent-based model for influenza transmission in Hong Kong
using the Covasim framework (https://github.com/InstituteforDiseaseModeling/covasim).
It simulates six influenza seasons (2009-2013) in Hong Kong and
evaluates various intervention strategies.

Key features:
- Age-structured population based on Hong Kong demographics
- Calibration to empirical data
- Multiple transmission layers (household, school, workplace, community)
- Vaccination interventions with varying efficacy
- Testing and self-isolation policies
- Mask-wearing interventions
- School holiday effects on transmission

Author: Liping Peng
Institution: The University of Hong Kong
Date: 9 September 2025

"""

import copy
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import covasim as cv
import os
import sciris as sc
import re
import base_func as bfc
from base_data import dis_prog
from base_subgroup import Sub_vac_cover, Sub_tip_base, Sub_tip_param




""" ************************* calibration *************************"""
def calib_change_beta(sim, calib_pars):
    change_beta = sim.get_intervention(cv.change_beta)
    change_beta.changes[0] = calib_pars['changes0']
    return sim


def main_calibrate(calib_interventions):
    # sim_calib = cv.Sim(pars_base, datafile=df_base, label="Baseline",
    #                    interventions=calib_interventions)
    with sc.timer('loading'):
        sim_calib = cv.Sim(pars_base, datafile=df_base, label="calibration",
                           interventions=calib_interventions,
                           popfile='result/load_people/my-people%s.ppl' % Season)


    calib_pars = dict(
        beta=[set_beta, set_beta-0.0015, set_beta+0.003],  # set_beta = 0.0073
        pop_infected=[set_infected, set_infected-2000, set_infected+10000],   # set_infected =  8000
        LP_sus_0_25  = [1, 0.5, 2.0],    # np.array([1.8, 1.8, 1.0, 1.0, 1.2])
        # LP_sus_25_45 = [1, 0.5, 1.5],
        LP_sus_45_65 = [1, 0.5, 2.0],
        LP_sus_65    = [1, 0.5, 2.0],
        LP_base_tip       = [0.05, 0.01, 0.10],

    )


    # 检查是否存在db,若存在则删除
    db_name = "result/result_main/calibration/%sLP_cali.db" % Season
    file_path = os.path.join(db_name)
    if os.path.exists(file_path):
        os.remove(file_path)

    calib = sim_calib.calibrate(calib_pars=calib_pars,
                                db_name=db_name,
                                # custom_fn=calib_change_beta,
                                # total_trials=20,
                                n_trials=n_trials,
                                n_workers=n_workers,
                                keep_db=True)



    print(calib.df)
    pl.ion()
    calib.plot_sims(to_plot=['cum_infections', 'new_infections'], do_save=True, fig_path="%s/figure/Figure_calib.pdf" % (output_path))
    if mismatch_age_bool:
        calib.LP_plot_age(is_save=True, save_path= "%s/figure/Figure_calib_byage.pdf" % (output_path))
    if mismatch_absent_bool:
        calib.LP_plot_absent(is_save=True, save_path="%s/figure/Figure_calib_absent.pdf" % (output_path))
    pl.pause(1)  ## show 1 second
    pl.close()  ## close fig, then continue running




""" ******************************Running simulation (output each sim's detail)******************************"""
def main_simulation(begin_seed, N):

    sims = []
    sims_seed = []
    r0_list = []
    re_list = []
    """ (1) single scenario and output model agent details"""
    for n in range(begin_seed, begin_seed + N):
        pars_base['beta'] = df_calib.loc[n,'beta']
        pars_base['pop_infected'] = int(df_calib.loc[n,'pop_infected'])
        pars_base['LP_sus_0_25'] = df_calib.loc[n, 'LP_sus_0_25']
        # pars_base['LP_sus_25_45'] = df_calib.loc[n, 'LP_sus_25_45']
        pars_base['LP_sus_45_65'] = df_calib.loc[n, 'LP_sus_45_65']
        pars_base['LP_sus_65'] = df_calib.loc[n, 'LP_sus_65']
        pars_base['LP_base_tip']  = df_calib.loc[n, 'LP_base_tip']


        # ## option 1: create population every time
        # sim_base = cv.Sim(pars_base, datafile=df_base, label=scenario,
        #                   # analyzers=[anlzr_daily_layer(),anlzr_daily_ages()],
        #                   interventions=scen_intervention,
        #                   rand_seed=1,
        #                   # variants = wild
        #                   )  # ,interventions=child_sus  #analyzers=cv.daily_stats()
        # sim_base.run()

        ## option 2: create at first, then load the same population
        if n == 0:
            with sc.timer('creating'):
                sim1 = cv.Sim(pars_base, datafile=df_base, label=scenario,
                              # analyzers=[anlzr_daily_ages(),anlzr_daily_layer()],
                              interventions=scen_intervention,
                              rand_seed=1).init_people()
            sim1.people.save('result/load_people/my-people%s.ppl' % Season)
        with sc.timer('loading'):
            sim_base = cv.Sim(pars_base, datafile=df_base, label=scenario,
                              # analyzers=[anlzr_daily_ages(),anlzr_daily_layer()],
                              # analyzers=[anlzr_daily_ages()], ## 暂不输出by layer的信息，节省时间。若恢复使用，需要对应修改def analyzer_output(n)
                              interventions=scen_intervention,
                              rand_seed=1, popfile='result/load_people/my-people%s.ppl' % Season).init_people()
        sim_base.run()


        """save sim for further analysis"""
        sim_base.save('%s/save_sim/save_sim%s.sim' % (output_path, n), keep_people=True)
        """ output details of each sim """
        sim_base.to_excel("%s/sims_detail/sim%s.xlsx" % (output_path, n))





if __name__ == "__main__":
    # Set working directory
    path = "D:/ABM influenza/"
    os.chdir(path)


    # =============================================================================
    # SIMULATION PARAMETERS
    # =============================================================================
    Season = 1

    # Vaccine parameters
    flu_vaccine = 'lp_flu_30'  # mismatch: "lp_flu_30"; match: "lp_flu_70"
    flu_vac_eff = float(re.search("lp_flu_(\d+)", flu_vaccine)[1]) * 0.01

    # Biological parameters
    imprt_bool = False  # Whether to include immune imprinting
    preHAI_bool = True  # Whether to include pre-existing immunity
    Flutype = 'H1N1'    # Influenza type: 'H1N1' for Season 1,3,5; 'H3N2' for Season 2,4,6


    # Calibration targets and simulation setting
    run_calibration = False  # False: run simulation, True: run calibration
    mismatch_age_bool = True  # Include age-specific infection data in calibration
    mismatch_absent_bool = True  # Include absenteeism data in calibration
    if run_calibration:
        n_trials = 8000
        n_workers = 5
    else:
        df_calib = pd.read_excel("result/result_main/calibration/Season%s_calib_param.xlsx" % (Season))  # Load calibrated parameters for simulation
        begin_seed = 0  # Starting random seed
        run_times = 1  # Number of simulation runs



    # Define key dates for the simulation timeline
    start_date = '2009-06-20'  # Model simulation start date
    end_date = '2009-11-21'  # Model simulation end date
    actual_start = '2009-07-05' # Actual outbreak start date
    actual_end = '2010-01-16'  # Actual outbreak end date
    press_date = '2009-06-27'  # Date of official press release (source: HK CHP)
    school_holiday = True  # Whether to account for school holiday effects
    school_year = '2009'   # Academic year for holiday scheduling


    # Hong Kong population structure by age group (here uses 2009 as an example)
    HK_pop = {
            '0-4': 226000,
            '5-9': 262900,
            '10-19':806000,
            '20-29':986600,
            '30-39':1107500,
            '40-49':1267700,
            '50-59':1082100,
            '60-69':559000,
            '70-79':435000,
            '80+':240000,}
    total_pop = sum([HK_pop[key] for key in HK_pop])
    # cv.data.country_age_data.data['China, Hong Kong Special Administrative Region'] = HK_pop   # Set Hong Kong population data for synthetic population generation. Only used with 'hybrid' or 'random' pop_type, not with 'synthpops'
    set_pop_scale = 100   # Population scaling factor (for computational efficiency)



    # =============================================================================
    # LOAD AND PROCESS EMPIRICAL DATA FOR CALIBRATION
    # =============================================================================
    ### Load observed epidemiological data for calibration target
    ### This data serves as the ground truth for model calibration
    df_base = pd.read_excel('LP_input data/data_season%s.xlsx' % Season, sheet_name="data")

    # Determine whether to include age-specific and absenteeism data in calibration
    if mismatch_age_bool:
        if mismatch_absent_bool:
            df_base = df_base[['date', 'new_infections', "new_infections.0.24", "new_infections.25.44", "new_infections.45.64", "new_infections...65", "new_absent"]]
        else:
            df_base = df_base[['date', 'new_infections', "new_infections.0.24", "new_infections.25.44", "new_infections.45.64","new_infections...65"]]
    else:
        df_base = df_base[['date', 'new_infections']]


    # Apply moving average to smooth observed data and reduce noise
    rolling_bool = True  # Enable/disable moving average smoothing
    roll_absent = True   # Whether to apply smoothing to absenteeism data
    if rolling_bool:
        if not roll_absent:
            cols_to_roll = df_base.columns.difference(['date', 'new_absent'])
            df_base[cols_to_roll] = df_base[cols_to_roll].rolling(window=7, center=True, min_periods=1).mean()
        else:
            df_base.iloc[:, 1:] = df_base.iloc[:, 1:].rolling(window=7, center=True, min_periods=1).mean() ## rolling avearge on 7 day window  ## 除了第一列date，都rolling mean
    print (df_base)

    # Restrict calibration data to the model simulation period only
    df_base = df_base[(df_base['date'] >= start_date) & (df_base['date'] <= end_date)]




    # =============================================================================
    # Calibrated parameters and optimal values
    # =============================================================================
    # Per-contact transmissibility (beta)
    set_beta = 0.0073
    # Initial infected population
    set_infected =  8000



    # =============================================================================
    # INTERVENTION SCENARIOS DEFINITION
    # =============================================================================
    ## Define intervention scenarios to simulate
    ## Scenario naming convention: S{scenario_id}_vg{vacc_group}_vc{vacc_coverage}_tip{prob_stay_home}_tie{stay_home_effect}_tida{early_termination}_mc{mask_coverage}_me{mask_effect}_intvd{interv_timing}_tidd{test_delay}
    ## Parameter details:
    ## - vg: Vaccination target group (refer to base_subgroup.py)
    ## - vc: Vaccine coverage increase (percentage points above baseline coverage)
    ## - tip: Daily probability of adherecing to stay-at-home among healthcare-seeking cases. (% above baseline)
    ## - tie: Stay-at-home effectiveness (%)
    ## - tida: Early termination of stay-at-home (days before recovery to end stay-at-home)
    ## - mc: Mask wearing coverage increase (% above baseline coverage)
    ## - me: Mask effectiveness (currently fixed at baseline)
    ## - intvd: Intervention timing (weeks after press release to implement, 0 = immediately)
    ## - tidd: Testing delay (days from symptom onset to diagnosis, 99 = no delay)
    ## Example: 'S0_vg0_vc0_tip0_tie0_tida0_mc0_me0_intvd0_tidd0' = Baseline scenario with no interventions

    list_Scenario = [
        'S0_vg0_vc0_tip0_tie0_tida0_mc0_me0_intvd0_tidd0',  # Baseline scenario (no additional interventions)
        # 'S1_vg1_vc10_tip0_tie0_tida0_mc0_me0_intvd0_tidd0'  # Example intervention scenario. vg1: Expand vaccination to age group 0-11 years; vc10: Increase vaccination coverage by 10 percentage points on top of baseline
    ]



    # =============================================================================
    # SCENARIO LOOP - EXECUTE SIMULATION FOR EACH INTERVENTION SCENARIO
    # =============================================================================
    for scenario in list_Scenario:
        output_path = 'result/result_main/Season%s/%s' % (Season, scenario)

        # Create output directory path for this scenario
        bfc.func_create_folder(output_path, 'figure')         # For output figures and plots
        bfc.func_create_folder(output_path, 'sims_detail')    # For detailed simulation results (Excel files)
        bfc.func_create_folder(output_path, 'sims_multisim')  # For multi-simulation ensemble results
        bfc.func_create_folder(output_path, 'save_sim')       # For saved simulation objects (.sim files)


        # =========================================================================
        # BASELINE INTERVENTION PARAMETERS
        # =========================================================================
        base_tie = 0.8       # Baseline stay-at-home effectiveness (80% reduction in transmission for stay-at-home individuals in non-household settings)
        base_tidd = 1        # Default testing delay = 1 day (time from symptom onset to diagnosis)
        base_mc = 0.10       # Baseline mask wearing coverage (10% of population wears masks)
        base_me = 0.5 * 0.5  # Baseline mask effectiveness: 50% efficacy × 50% compliance = 25% overall effect



        # =========================================================================
        # PARSE INTERVENTION PARAMETERS FROM SCENARIO STRING
        # =========================================================================
        # Extract and convert parameter values from matched groups
        param_value = re.search("S\d+_(vg\d+)_vc(\d+)_tip(\d+)_tie(\d+)_tida(\d+)_mc(\d+)_me(\d+)_intvd(\d+)_tidd(\d+)", scenario)
        param_vg    = param_value.group(1)
        param_vc    = float(param_value.group(2)) * 0.01
        param_tip   = float(param_value.group(3)) * 0.01
        param_tida  = int(param_value.group(5))
        # use baseline if not specified in scenario
        if float(param_value.group(4)) > 0:
            param_tie  = float(param_value.group(4)) * 0.01
        else:
            param_tie  = base_tie
        # Testing delay handling - special case: 99 means no delay (delay = 0)
        if int(param_value.group(9)) >0:
            if int(param_value.group(9)) == 99:   # Special code: 99 = no testing delay
                param_tidd = 0
            else:
                param_tidd  = int(param_value.group(9))
        else:
            param_tidd  = base_tidd   # default delay = 1day
        param_mc    = base_mc + float(param_value.group(6)) * 0.01
        param_me    = base_me
        # Intervention timing - when to implement the intervention package
        param_intvd = int(param_value.group(8))   # Weeks after press release to implement
        if param_intvd == 0:
            intv_date = press_date   # Implement on press release date
        elif param_intvd >= 1:
            # Calculate implementation date: n weeks after press release
            intv_date = datetime.datetime.strptime(press_date, "%Y-%m-%d") + datetime.timedelta(days= param_intvd * 7)
            intv_date = intv_date.strftime("%Y-%m-%d")




        # =========================================================================
        # SEASONAL TRANSMISSION PATTERN (SINUSOIDAL MODEL)
        # =========================================================================
        ### This must be set before mask interventions as both modify beta values
        # Calculate simulation duration in days
        duration = datetime.datetime.strptime(end_date, "%Y-%m-%d") - datetime.datetime.strptime(start_date, "%Y-%m-%d")
        total_days = np.arange(duration.days)
        # Initialize beta values with baseline transmission (1.0 = 100%)
        beta_vals = [1] * duration.days
        # Calculate start date for seasonal pattern (July 15th)
        turndate1 = datetime.datetime.strptime("2009-07-15", "%Y-%m-%d") - datetime.datetime.strptime(start_date,
                                                                                                      "%Y-%m-%d")
        # Seasonal pattern parameters
        freq = 163
        amplit = 0.38
        phase = -17
        # Apply cosine function to create seasonal pattern after July 15th
        beta_vals[turndate1.days:] = (np.cos(
            2 * np.pi * (np.arange(duration.days - turndate1.days) + phase) / freq) * amplit + 1).tolist()
        beta_vals = np.array(beta_vals)




        # =========================================================================
        # DEFINE INTERVENTION FUNCTIONS
        # =========================================================================
        ### (1) INTERVENTION: Vaccination
        intv_vac = cv.vaccinate_prob(vaccine=flu_vaccine,
                                     days=1,                  # Start on day 1 (pre-season vaccination campaign)
                                     subtarget=Sub_vac_cover  # Use custom targeting function for age-specific coverage
                                     )


        ### (2) INTERVENTION: Stay-at-home
        if run_calibration:
            # For calibration: apply throughout simulation
            base_TI = cv.test_prob(
                symp_prob=1,      # 100% testing probability for symptomatic individuals
                asymp_prob=0,     # 0% testing probability for asymptomatic individuals
                start_day=1,      # Start testing from day 1 of simulation
                test_delay=param_tidd,  # Days from symptom onset to diagnosis
                subtarget=Sub_tip_base  # Custom function for baseline daily probabilities
            )
        elif param_tip <= 0 :
            # For scenarios with no stay-at-home improvement: use baseline setting
            base_TI = cv.test_prob(symp_prob=1,
                                   asymp_prob=0,
                                   start_day=1,
                                   test_delay=param_tidd,
                                   subtarget=Sub_tip_base)
        else:
            # For intervention scenarios: two-phase strategy
            # Phase 1: Baseline setting until intervention start date
            iso_start = datetime.datetime.strptime(intv_date, "%Y-%m-%d") - datetime.datetime.strptime(start_date, "%Y-%m-%d")
            base_TI = cv.test_prob(symp_prob=1,
                                   asymp_prob=0,
                                   start_day=1,
                                   end_day= iso_start.days-1,
                                   test_delay=param_tidd,
                                   subtarget=Sub_tip_base)
            # Phase 2: Enhanced setting after intervention start date
            intv_TI = cv.test_prob(symp_prob=1,
                                   asymp_prob=0,
                                   start_day=iso_start.days,
                                   test_delay=param_tidd,
                                   subtarget=Sub_tip_param)



        ### (3) INTERVENTION: Mask
        # Calculate mask effectiveness as an odds ratio
        mask_OR = 1 - param_me * param_mc

        # Apply mask intervention to s,w,c layers  after start date
        mask_start = datetime.datetime.strptime(intv_date, "%Y-%m-%d") - datetime.datetime.strptime(start_date,"%Y-%m-%d")
        temp = [1] * len(total_days)
        temp[mask_start.days:duration.days] = [mask_OR] * (duration.days - mask_start.days)
        beta_vals_mask = beta_vals * temp

        ### (3.1): change beta for household (seasonal pattern only)
        beta_h = copy.deepcopy(beta_vals)
        intv_cb_h = cv.change_beta(days=total_days, changes=beta_h, do_plot=False, layers=['h'])  ## beta_layer  = dict(h=3.0, s=0.6, w=0.6, c=0.3), 已经设置了不同layer的beta weight， 所以在实施mask的时候，不区分layer？


        ### (3.2): change beta for school (seasonal pattern + mask + school holiday)
        beta_s = copy.deepcopy(beta_vals_mask)
        if school_holiday:
            # Account for school holiday effects
            temp_holi = [1] * len(total_days)
            pre_holi_start = datetime.datetime.strptime(start_date, "%Y-%m-%d") - datetime.datetime.strptime(start_date, "%Y-%m-%d")  # Partial closure before holidays (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3206396/)
            holi_start = datetime.datetime.strptime('%s-07-15' % school_year, "%Y-%m-%d") - datetime.datetime.strptime(start_date, "%Y-%m-%d")
            holi_end = datetime.datetime.strptime('%s-08-31' % school_year, "%Y-%m-%d") - datetime.datetime.strptime(start_date, "%Y-%m-%d")
            # Only holiday reduction (no mask effect during holidays)
            beta_s[pre_holi_start.days:holi_start.days] = [0.5] * (holi_start.days - pre_holi_start.days) # Partial closure before holidays
            beta_s[holi_start.days:holi_end.days] = [0.2] * (holi_end.days - holi_start.days)             # Full closure during holidays
        intv_cb_s = cv.change_beta(days=total_days, changes=beta_s, do_plot=False, layers=['s'])  # Apply to school layer


        ### (3.3): change beta for workplace and community (seasonal pattern + mask)
        beta_wc = copy.deepcopy(beta_vals_mask)
        intv_cb_wc = cv.change_beta(days=total_days, changes=beta_wc, do_plot=False, layers=['w', 'c'])


        ### COMBINE INTERVENTIONS INTO FINAL SET
        if run_calibration:
            scen_intervention = [intv_vac, base_TI,          intv_cb_h, intv_cb_s, intv_cb_wc]
        elif param_tip <= 0:
            scen_intervention = [intv_vac, base_TI,          intv_cb_h, intv_cb_s, intv_cb_wc]
        else:
            scen_intervention = [intv_vac, base_TI, intv_TI, intv_cb_h, intv_cb_s,intv_cb_wc]



        # =========================================================================
        # PARAMETERS SET FOR SIMULATIONS
        # =========================================================================
        pars_base = dict(
            location='Hong Kong',  # Location
            pop_type='synthpops',  # Population generation method: synthpops, hybrid, or random
            LP_season=Season,      # Record Season
            LP_flutype=Flutype,    # Record variant
            LP_imprinting=imprt_bool,  # Whether to include immune imprinting effects
            LP_preHAI=preHAI_bool,     # Whether to include pre-existing immunity

            # Calibration targets
            LP_mismatch_age=mismatch_age_bool,  # Include age-specific data in calibration
            LP_mismatch_absent=mismatch_absent_bool,  # Include absenteeism data in calibration

            # Calibrated parameters (will be overwritten during calibration)
            beta=set_beta,  # beta
            LP_sus_0_25=1,  # Susceptibility for 0-25 age group
            LP_sus_25_45=1,  # Baseline susceptibility group (other groups relative to this)
            LP_sus_45_65=1,  # Susceptibility for 45-65 age group
            LP_sus_65=1,  # Susceptibility for 65+ age group
            LP_base_tip=0,  # Baseline daily probability of adherence to stay at home

            # Population parameters
            pop_scale=set_pop_scale,  # Scale factor (100 = 1% of real population)
            pop_size=round(total_pop / set_pop_scale, 0),  # Scaled population size

            # Simulation timing
            start_day=start_date,  # Simulation start date
            end_day=end_date,  # Simulation end date

            # transmission
            pop_infected=set_infected,  # Initial number of infected individuals
            dur=dis_prog,  # Disease progression parameters
            LP_symp_probs=np.array([0.80, 0.80, 0.80, 0.80, 0.80]),  # Symptomatic probabilities by age group
            LP_trans_ORs=np.array([2.00, 1.00, 1.00, 1.00, 1.00], dtype=float),  # Transmission odds ratios by age

            # Intervention
            LP_param_vg=param_vg,  # Vaccination target group
            LP_param_vc=param_vc,  # Vaccination coverage increase
            LP_param_tip=param_tip,  # Stay-at-home probability increase
            iso_factor=dict(h=1.0, s=1 - param_tie, w=1 - param_tie, c=1 - param_tie, l=1.0),  # Stay-at-home effectiveness per layer
            LP_end_tida=param_tida,  # Early termination of stay-at-home (days before recovery)

        )



        # =========================================================================
        # EXECUTE MAIN SIMULATION PROCESS
        # =========================================================================
        if 'S0_' in scenario:  # Baseline scenario
            if run_calibration:
                main_calibrate(calib_interventions=scen_intervention)   # Run calibration
            else:
                main_simulation(begin_seed, run_times)                  # Run simulation with calibrated parameters
        else:
            main_simulation(begin_seed, run_times)                      # Run simulation with intervention parameters






