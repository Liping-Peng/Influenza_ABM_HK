"""
This module contains data for all seasons.

Author: Liping Peng
Institution: The University of Hong Kong
Date: 9 September 2025
"""


# =============================================================================
# COMMON DATA FOR ALL SEASONS
# =============================================================================

### Disease Progression Parameters
dis_prog = dict(exp2inf=dict(dist='lognormal_int', par1=1.9, par2=0.2),
                # Ref: https://www.thelancet.com/journals/lancet/article/PIIS1473-3099(09)70069-6/fulltext
                # Ref: https://www.nejm.org/doi/10.1056/NEJMoa0906089?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200www.ncbi.nlm.nih.gov

                # Infectiousness to symptom onset (negligible for influenza)
                inf2sym=dict(dist='lognormal_int', par1=0, par2=0),

                # Assumed same duration as symptomatic cases
                asym2rec=dict(dist='lognormal_int', par1=6.0, par2=1.02),
                mild2rec=dict(dist='lognormal_int', par1=6.0, par2=1.02),

                # Currently severe/critical probabilities = 0, so this duration is less critical
                # Ref: https://www.nejm.org/doi/10.1056/NEJMoa0906089?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200www.ncbi.nlm.nih.gov
                sym2sev=dict(dist='lognormal_int', par1=6.6, par2=4.9),
                sev2crit=dict(dist='lognormal_int', par1=1.5, par2=2.0),
                sev2rec=dict(dist='lognormal_int', par1=18.1, par2=6.3),
                crit2rec=dict(dist='lognormal_int', par1=18.1, par2=6.3),
                crit2die=dict(dist='lognormal_int', par1=10.7, par2=4.8),
                )




dict_ModelStart = {1: '2009-06-20', 2: '2010-06-25', 3: '2010-12-04', 4: '2012-03-01', 5: '2013-01-01', 6: '2013-05-31', }
dict_ModelEnd = {1: '2009-11-21', 2: '2011-01-16', 3: '2011-05-01', 4: '2012-09-23', 5: '2013-07-20', 6: '2013-12-12', }

dict_ActualStart = {1:'2009-07-05', 2:"2010-08-01", 3:"2011-01-09", 4:"2012-03-11", 5:"2013-02-17", 6:"2013-07-14",}
dict_ActualEnd = {1:"2009-11-21", 2:"2010-10-16", 3:"2011-02-26", 4:"2012-06-23", 5:"2013-04-20", 6:"2013-10-12",}

dict_press = {1:'2009-06-27', 2:'2010-07-29',3:'2011-01-17',4:'2012-03-26',5:'2013-02-04',6:'2013-07-22'}

dict_pop = {1:6972800, 2:7024200, 3:7052100, 4:7150100, 5:7171000, 6:7178900}

dict_pop_byage = {1: {0:1788200,	25:2234650,	  45:1995450,	  65:954500,},
                  2: {0:1765500,	25:2225950,   45:2040100,	  65:992650,},
                  3: {0:1750950,	25:2219050,	  45:2065500,	  65:1016600,},
                  4: {0:1725900,	25:2227150,	  45:2135100,	  65:1061950,},
                  5: {0:1709850,	25:2217150,	  45:2162250,	  65:1081750,},
                  6: {0:1689200,	25:2211250,	  45:2184800,	  65:1093650,}
                  }


dict_vg_label = {0: 'Baseline',
                 1: 'G1: <12',
                 2: 'G2: <18',
                 3: 'G3: >=65',
                 4: 'G4: <18 or >=65',
                 5: 'G5: universal'
                 }

list_vg_label = ['G1: <12', 'G2: <18', 'G3: >=65', 'G4: <18 or >=65']



dict_school_closure = {1: {'start':'2009-07-15', 'end':'2009-08-31'},
                       2: {'start':'2010-07-15', 'end':'2010-08-31'},
                       3: {'start':'2011-01-10', 'end':'2011-01-23'},
                       4: {'start':'2012-03-11', 'end':'2012-03-24'},
                       5: {'start':'2013-02-20', 'end':'2013-03-05'},
                       6: {'start':'2013-07-15', 'end':'2013-08-31'}
                       }

