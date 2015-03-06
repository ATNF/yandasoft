#! /usr/bin/env python

from antennas import *

datafile = "ConfigurationData/ASKAP-SEIC-0005_Antenna_Configuration.csv"

def getConfigs():

    config = {}
    config["A27CR3P6B"]                = AntennaConfig("A27CR3P6B",                range(1, 37),        datafile)
    config["ASKAP36"]                  = AntennaConfig("ASKAP36",                  range(1, 37),        datafile)
    config["A27CR3"]                   = AntennaConfig("A27CR3",                   range(1, 31),        datafile)
    config["BETA"]                     = AntennaConfig("BETA",                     (1, 3, 6, 8, 9, 15), datafile)
    config["BETAXYZ"]                  = AntennaConfig("BETAXYZ",                  (1, 3, 6, 8, 9, 15), datafile)
    config["BETA2"]                    = AntennaConfig("BETA2",                    (1, 3, 6, 8, 9,  2), datafile)
    config["BETA4"]                    = AntennaConfig("BETA4",                    (1, 3, 6, 8, 9,  4), datafile)
    config["BETA5"]                    = AntennaConfig("BETA5",                    (1, 3, 6, 8, 9,  5), datafile)
    config["BETA7"]                    = AntennaConfig("BETA7",                    (1, 3, 6, 8, 9,  7), datafile)
    config["BETA10"]                   = AntennaConfig("BETA10",                   (1, 3, 6, 8, 9, 10), datafile)
    config["BETA11"]                   = AntennaConfig("BETA11",                   (1, 3, 6, 8, 9, 11), datafile)
    config["BETA12"]                   = AntennaConfig("BETA12",                   (1, 3, 6, 8, 9, 12), datafile)
    config["BETA13"]                   = AntennaConfig("BETA13",                   (1, 3, 6, 8, 9, 13), datafile)
    config["BETA14"]                   = AntennaConfig("BETA14",                   (1, 3, 6, 8, 9, 14), datafile)
    config["BETA15"]                   = AntennaConfig("BETA15",                   (1, 3, 6, 8, 9, 15), datafile)
    config["BETA16"]                   = AntennaConfig("BETA16",                   (1, 3, 6, 8, 9, 16), datafile)
    config["BETA17"]                   = AntennaConfig("BETA17",                   (1, 3, 6, 8, 9, 17), datafile)
    config["BETA18"]                   = AntennaConfig("BETA18",                   (1, 3, 6, 8, 9, 18), datafile)
    config["BETA19"]                   = AntennaConfig("BETA19",                   (1, 3, 6, 8, 9, 19), datafile)
    config["BETA20"]                   = AntennaConfig("BETA20",                   (1, 3, 6, 8, 9, 20), datafile)
    config["BETA21"]                   = AntennaConfig("BETA21",                   (1, 3, 6, 8, 9, 21), datafile)
    config["BETA22"]                   = AntennaConfig("BETA22",                   (1, 3, 6, 8, 9, 22), datafile)
    config["BETA23"]                   = AntennaConfig("BETA23",                   (1, 3, 6, 8, 9, 23), datafile)
    config["BETA24"]                   = AntennaConfig("BETA24",                   (1, 3, 6, 8, 9, 24), datafile)
    config["BETA25"]                   = AntennaConfig("BETA25",                   (1, 3, 6, 8, 9, 25), datafile)
    config["BETA26"]                   = AntennaConfig("BETA26",                   (1, 3, 6, 8, 9, 26), datafile)
    config["BETA27"]                   = AntennaConfig("BETA27",                   (1, 3, 6, 8, 9, 27), datafile)
    config["BETA28"]                   = AntennaConfig("BETA28",                   (1, 3, 6, 8, 9, 28), datafile)
    config["BETA29"]                   = AntennaConfig("BETA29",                   (1, 3, 6, 8, 9, 29), datafile)
    config["BETA29A"]                  = AntennaConfig("BETA29A",                  (1, 3, 6, 8,11, 29), datafile)
    config["BETA29B"]                  = AntennaConfig("BETA29B",                  (1, 3, 8, 9,11, 29), datafile)
    config["BETA30"]                   = AntennaConfig("BETA30",                   (1, 3, 6, 8, 9, 30), datafile)
    config["BETA31"]                   = AntennaConfig("BETA31",                   (1, 3, 6, 8, 9, 31), datafile)
    config["BETA32"]                   = AntennaConfig("BETA32",                   (1, 3, 6, 8, 9, 32), datafile)
    config["BETA33"]                   = AntennaConfig("BETA33",                   (1, 3, 6, 8, 9, 33), datafile)
    config["BETA34"]                   = AntennaConfig("BETA34",                   (1, 3, 6, 8, 9, 34), datafile)
    config["BETA35"]                   = AntennaConfig("BETA35",                   (1, 3, 6, 8, 9, 35), datafile)
    config["BETA36"]                   = AntennaConfig("BETA36",                   (1, 3, 6, 8, 9, 36), datafile)
    config["ADE6_POSS_SJ1"]            = AntennaConfig("ADE6_POSS_SJ1",            ( 4, 5,11,17,25,26), datafile)
    config["ADE6_POSS_SJ2"]            = AntennaConfig("ADE6_POSS_SJ2",            ( 4, 5,11,12,19,20), datafile)
    config["ADE6_POSS_MJK"]            = AntennaConfig("ADE6_POSS_MJK",            ( 4, 5,10,11,18,19), datafile)
    config["ADE6_POSS_JB"]             = AntennaConfig("ADE6_POSS_JB",             (10,11,18,21,22,29), datafile)
    config["ADE6_POSS_JB2"]            = AntennaConfig("ADE6_POSS_JB",             (10,11,18,21,22, 7), datafile)
    config["ADE6_POSS_JB3"]            = AntennaConfig("ADE6_POSS_JB",             (10,11,18,21,22,19), datafile)
    config["ADE12_POSS_EMU1.5"]        = AntennaConfig("ADE12_POSS_EMU1.5",        ( 4, 5, 7,10,11,12,14,16,18,21,23,26), datafile)
    config["ADE12_POSS_EMU3"]          = AntennaConfig("ADE12_POSS_EMU3",          ( 4, 5,10,18,19,20,25,27,28,29,35,36), datafile)
    config["ADE12_POSS_EMU6"]          = AntennaConfig("ADE12_POSS_EMU6",          ( 4, 5,10,17,19,20,27,29,32,33,35,36), datafile)
    config["ADE12_POSS_WALLABY"]       = AntennaConfig("ADE12_POSS_WALLABY",       ( 2, 4, 5, 7,10,11,12,14,16,17,18,19), datafile)
    config["ADE6_AK07_PROPOSED"]       = AntennaConfig("ADE6_AK07_PROPOSED",       ( 2, 4, 5, 7,10,12), datafile)
    config["ADE6_AK19_PROPOSED"]       = AntennaConfig("ADE6_AK19_PROPOSED",       ( 2, 4, 5,19,10,12), datafile)
    config["ADE6_AK11_PROPOSED"]       = AntennaConfig("ADE6_AK11_PROPOSED",       ( 2, 4, 5,11,10,12), datafile)
    config["ADE6_AK14_PROPOSED"]       = AntennaConfig("ADE6_AK14_PROPOSED",       ( 2, 4, 5,14,10,12), datafile)
    config["ADE6_AK17_PROPOSED"]       = AntennaConfig("ADE6_AK17_PROPOSED",       ( 2, 4, 5,17,10,12), datafile)
    config["ADE12_AK07_PROPOSED"]      = AntennaConfig("ADE12_AK07_PROPOSED",      ( 2, 4, 5, 7,10,12,24,27,30,13,16,28), datafile)
    config["ADE12_AK19_PROPOSED"]      = AntennaConfig("ADE12_AK19_PROPOSED",      ( 2, 4, 5,19,10,12,24,27,30,13,16,28), datafile)
    config["ADE12_AK11_PROPOSED"]      = AntennaConfig("ADE12_AK11_PROPOSED",      ( 2, 4, 5,11,10,12,24,27,30,13,16,28), datafile)
    config["ADE12_AK14_PROPOSED"]      = AntennaConfig("ADE12_AK14_PROPOSED",      ( 2, 4, 5,14,10,12,24,27,30,13,16,28), datafile)
    config["ADE12_AK17_PROPOSED"]      = AntennaConfig("ADE12_AK17_PROPOSED",      ( 2, 4, 5,17,10,12,24,27,30,13,16,28), datafile)
    config["ADE18_AK07_PROPOSED"]      = AntennaConfig("ADE18_AK07_PROPOSED",      ( 2, 4, 5, 7,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
    config["ADE18_AK11_PROPOSED"]      = AntennaConfig("ADE18_AK11_PROPOSED",      ( 2, 4, 5,11,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
    config["ADE18_AK21_PROPOSED"]      = AntennaConfig("ADE18_AK21_PROPOSED",      ( 2, 4, 5,21,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
    config["ADE18_AK23_PROPOSED"]      = AntennaConfig("ADE18_AK23_PROPOSED",      ( 2, 4, 5,23,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
    config["ADE18_AK17_PROPOSED"]      = AntennaConfig("ADE18_AK17_PROPOSED",      ( 2, 4, 5,17,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
    config["ADE12_WALLABY_A"]          = AntennaConfig("ADE12_WALLABY_A",          ( 1, 4, 5, 6, 7,11,13,14,15,19,21,25), datafile)
    config["ADE12_WALLABY_B"]          = AntennaConfig("ADE12_WALLABY_B",          ( 1, 2, 6, 9,12,15,16,18,23,27,28,29), datafile)
    config["ADE12_WALLABY_C"]          = AntennaConfig("ADE12_WALLABY_C",          ( 2, 3, 4, 6, 7,12,15,18,19,22,23,25), datafile)
    #
    # These are the configurations decided upon in 2013 for the 6-, 12- and 18-antenna configurations for the MkII PAFs
    config["ADE6"]                     = AntennaConfig("ADE6",                     ( 2, 4, 5,14,10,12), datafile)
    config["ADE12"]                    = AntennaConfig("ADE12",                    ( 2, 4, 5,14,10,12,24,27,30,13,16,28), datafile)
    config["ASKAP12"]                  = AntennaConfig("ASKAP12",                  ( 2, 4, 5,14,10,12,24,27,30,13,16,28), datafile)
    config["ADE18"]                    = AntennaConfig("ADE18",                    ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26), datafile)
    config["ASKAP18"]                  = AntennaConfig("ASKAP18",                  ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26), datafile)
    #
    # An alternative ADE12 config with AK34 replacing AK30 (to test the long fibre run)
    config["ADE12_AK34"]               = AntennaConfig("ADE12_AK34",               ( 2, 4, 5,14,10,12,24,27,34,13,16,28), datafile)
    #
    # The following are those under consideration for the ASKAP community meeting, October 21, 2014
    # 16-ant config with Batch #4 from Adrian Rispler's talk at October 21, 2014 community meeting - was ADE12_BATCH4
    config["ASKAP16_A"]                = AntennaConfig("ASKAP16_A",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19), datafile)
    # 16-ant config with Batch #5 from Adrian Rispler's talk at October 21, 2014 community meeting - was ADE12_BATCH5
    config["ASKAP16_B"]                = AntennaConfig("ASKAP16_B",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,23,26,11,17), datafile)
    # 16-antenna configs (ADE12 + 4, either 4BETA or 4 nonBETA)
    config["ADE12_4BETA_1368"]         = AntennaConfig("ADE12_4BETA_1368",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6, 8), datafile)
    config["ADE12_4BETA_1369"]         = AntennaConfig("ADE12_4BETA_1369",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6, 9), datafile)
    config["ADE12_4BETA_136F"]         = AntennaConfig("ADE12_4BETA_136F",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6,15), datafile)
    config["ADE12_4BETA_139F"]         = AntennaConfig("ADE12_4BETA_139F",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 9,15), datafile)
    config["ADE12_4NONBETA_1"]         = AntennaConfig("ADE12_4NONBETA_1",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,19,23,26,21), datafile)
    config["ADE12_4NONBETA_2"]         = AntennaConfig("ADE12_4NONBETA_2",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,17,18,20,23), datafile)
    config["ADE12_4NONBETA_3"]         = AntennaConfig("ADE12_4NONBETA_3",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,17,11,22,23), datafile)
    config["ADE12_4NONBETA_4"]         = AntennaConfig("ADE12_4NONBETA_4",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,26,11,22, 7), datafile)
    config["ADE12_4NONBETA_5"]         = AntennaConfig("ADE12_4NONBETA_5",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,19,23,26,11), datafile)
    config["ADE12_4NONBETA_6"]         = AntennaConfig("ADE12_4NONBETA_6",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,17), datafile)
    # These 16-antenna configs have a mix of BETA & non-BETA                       
    config["ADE12_2BETA_13_AK19_23"]   = AntennaConfig("ADE12_2BETA_13_AK19_23",   ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3,19,23), datafile)
    config["ASKAP16_C"]                = AntennaConfig("ASKAP16_C",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3,19,23), datafile)
    config["ADE12_2BETA_16_AK19_23"]   = AntennaConfig("ADE12_2BETA_16_AK19_23",   ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 6,19,23), datafile)
    config["ASKAP16_D"]                = AntennaConfig("ASKAP16_D",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 6,19,23), datafile)
    config["ADE12_2BETA_16_AK26_23"]   = AntennaConfig("ADE12_2BETA_16_AK26_23",   ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 6,26,23), datafile)
    config["ASKAP16_E"]                = AntennaConfig("ASKAP16_E",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 6,26,23), datafile)
    config["ADE12_2BETA_1F_AK19_23"]   = AntennaConfig("ADE12_2BETA_1F_AK19_23",   ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1,15,19,23), datafile)
    config["ASKAP16_F"]                = AntennaConfig("ASKAP16_F",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1,15,19,23), datafile)
    config["ADE12_1BETA_1_AK19_23_26"] = AntennaConfig("ADE12_1BETA_1_AK19_23_26", ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1,19,23,26), datafile)

    # 20-ant array using batches 4 & 5 - was ADE12_BATCH4_5
    config['ASKAP20_A']                = AntennaConfig("ASKAP20_A",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,11,17), datafile)
    # 20-antenna configs (ADE12 + 8, either 6BETA+2, or 8 non-BETA)
    config["ADE12_6BETA_AK31_33"]      = AntennaConfig("ADE12_6BETA_AK31_33",      ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,31,33), datafile)
    config["ADE12_6BETA_AK19_23"]      = AntennaConfig("ADE12_6BETA_AK19_23",      ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,19,23), datafile)
    config["ADE12_6BETA_AK26_33"]      = AntennaConfig("ADE12_6BETA_AK26_33",      ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,26,33), datafile)
    config["ADE12_6BETA_AK31_23"]      = AntennaConfig("ADE12_6BETA_AK31_23",      ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,31,23), datafile)
    config["ADE12_8NONBETA_1"]         = AntennaConfig("ADE12_8NONBETA_1",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,22,34), datafile)
    config["ADE12_8NONBETA_2"]         = AntennaConfig("ADE12_8NONBETA_2",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,17,21), datafile)
    config["ADE12_8NONBETA_3"]         = AntennaConfig("ADE12_8NONBETA_3",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,32,36), datafile)
    config["ADE12_8NONBETA_4"]         = AntennaConfig("ADE12_8NONBETA_4",         ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,11,17), datafile)
    # These 20-antenna configs have ADE-18 plus two BETA antennas:                 
    config["ADE18_2BETA_13"]           = AntennaConfig("ADE18_2BETA_13",           ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1, 3), datafile)
    config["ASKAP20_B"]                = AntennaConfig("ASKAP20_B",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1, 3), datafile)
    config["ADE18_2BETA_16"]           = AntennaConfig("ADE18_2BETA_16",           ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1, 6), datafile)
    config["ASKAP20_C"]                = AntennaConfig("ASKAP20_C",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1, 6), datafile)
    config["ADE18_2BETA_1F"]           = AntennaConfig("ADE18_2BETA_1F",           ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1,15), datafile)
    config["ASKAP20_D"]                = AntennaConfig("ASKAP20_D",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 1,15), datafile)
    config["ASKAP20_E"]                = AntennaConfig("ASKAP20_E",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 3, 6), datafile)
    config["ASKAP20_F"]                = AntennaConfig("ASKAP20_F",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 3, 9), datafile)
    config["ASKAP20_G"]                = AntennaConfig("ASKAP20_G",                ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26, 6, 9), datafile)

    # 30-antenna cases:
    a36=np.arange(1,37)
    # All 30 non-BETA antennas - was ASKAP30_NOBETA
    a=a36[(a36!=1)&(a36!=3)&(a36!=6)&(a36!=8)&(a36!=9)&(a36!=15)]
    config['ASKAP30_A']=AntennaConfig('ASKAP30_A',a.tolist(), datafile)
    # Put all the BETA ants in, and leave out current last two batches - was ASKAP30_BETA_BYBATCH
    a=a36[(a36!=25)&(a36!=7)&(a36!=18)&(a36!=32)&(a36!=36)&(a36!=29)]
    config['ASKAP30_B']=AntennaConfig('ASKAP30_B',a, datafile)
    # Replace 8,9,15 with 7,25,36 - was ASKAP30_3BETA
    a=a36[(a36!=15)&(a36!=8)&(a36!=18)&(a36!=32)&(a36!=9)&(a36!=29)]
    config['ASKAP30_C']=AntennaConfig('ASKAP30_C',a, datafile)
    # Replace 3,8,9,15 with 7,25,32,36 - was ASKAP30_2BETA
    a=a36[(a36!=15)&(a36!=8)&(a36!=18)&(a36!=3)&(a36!=9)&(a36!=29)]
    config['ASKAP30_D']=AntennaConfig('ASKAP30_D',a, datafile)
    # Replace 6,8,9,15 with 7,25,32,36
    a=a36[(a36!=15)&(a36!=8)&(a36!=18)&(a36!=6)&(a36!=9)&(a36!=29)]
    config['ASKAP30_E']=AntennaConfig('ASKAP30_E',a, datafile)
    # Replace 1,8,9,15 with 7,25,32,36
    a=a36[(a36!=15)&(a36!=8)&(a36!=18)&(a36!=1)&(a36!=9)&(a36!=29)]
    config['ASKAP30_F']=AntennaConfig('ASKAP30_F',a, datafile)
    # Alternative 3-BETA array: include 32 at expense of 11
    a=a36[(a36!=15)&(a36!=8)&(a36!=18)&(a36!=11)&(a36!=9)&(a36!=29)]
    config['ASKAP30_G']=AntennaConfig('ASKAP30_G',a, datafile)

    # These are 18-antenna configs suggested by WALLABY & DINGO at the community meeting
    config["ASKAP18_WALLABY"]          = AntennaConfig("ASKAP18_WALLABY",          ( 2, 3, 4, 5, 9,10,11,12,13,14,16,21,22,24,25,27,28,30), datafile)
    config["ASKAP18_DINGO"]            = AntennaConfig("ASKAP18_DINGO",            ( 2, 3, 4, 5, 6,10,11,12,13,14,15,16,21,24,25,27,28,30), datafile)

    return config

