#!/usr/bin/env python2.7

import csv
import constants
#import constants_one_energy_bin
import bayes
import numpy


OUT_DIR = "data_pass5f"
PATHS = \
        {"Segue"    :   {"EvList": "Ben/Pass5f/segue1_eventList_pass5f_wZnCorr.txt",  "Area": "Ben/Pass5d/segue1_eventList_pass5d_wZnCorr.txt.mce" }, \
         "Bootes1"  :   {"EvList": "Ben/Pass5f/bootes_eventList_pass5f_wZnCorr.txt",  "Area": "Ben/Pass5d/bootes_eventList_pass5d_wZnCorr.txt.mce"}, \
         "Draco"    :   {"EvList": "Ben/Pass5f/draco_eventList_pass5f_wZnCorr.txt",   "Area": "Ben/Pass5d/draco_eventList_pass5d_wZnCorr.txt.mce" }, \
         "UrsaMinor":   {"EvList": "Ben/Pass5f/umi_eventList_pass5f_wZnCorr.txt",     "Area": "Ben/Pass5d/umi_eventList_pass5d_wZnCorr.txt.mce"   }, \
         "WilmanI"  :   {"EvList": "Ben/Pass5f/willman_eventList_pass5f_wZnCorr.txt", "Area": "Ben/Pass5d/wilman_eventList_pass5d_wZnCorr.txt.mce"}}

Kfudge = 1.0

def GetAEff(filename):
    """ Returns RunInfo which is a dictionary with the following structure:
    {run1: {'E': [E1, E2, E3], 'A', [A1, A2, A3]}, run2: ...} """

    RunInfo = dict()
    Energies = []
    Areas = []

    def UpdateRunInfo(Dict, RunNumber, Energies, Areas):
        Dict[RunNumber] = {'E': Energies, 'A': Areas}

    PrevRun = -1
    RunCnt = 0
    for data in constants.get_event_line(filename, ' '):
        Run = int(data[0])
        if Run != PrevRun:
            if PrevRun != -1:
                UpdateRunInfo(RunInfo, PrevRun, Energies, Areas)
                Energies = []
                Areas = []
            RunCnt = RunCnt + 1
            PrevRun = Run
        Energy = data[2]
        Area = data[3]
        Energies.append(Energy)
        Areas.append(Area)
    UpdateRunInfo(RunInfo, PrevRun, Energies, Areas)
    return RunInfo


#for k, v in AreasInfo.iteritems():
#    print k
#    for i in range(0, len(v['E'])):
#        print i, v['E'][i], v['A'][i]
#    print "----------------------------------------------"

def Area(Energy, Runs, AreasInfo):
    return map(lambda Run: numpy.interp(Energy, AreasInfo[Run]['E'], AreasInfo[Run]['A']), Runs)

def Times(C, CRunArg, RunTimeArg):
    return map(lambda x: RunTimeArg[x], list(CRunArg[C]))

def GetRunsCAndTime(EventlistPath):
    """ Return RunC, CRun and RunTime (in same order) """
    RunC = dict()
    CRun = dict()
    RunTime = dict()
    # print "Assigning (pairing) weights to run numbers..."
    for data in constants.get_event_line(EventlistPath, ' '):
        # print data, data[9]
        Run = int(data[0])
        Alpha = data[9]
        if Alpha == "0": Alpha = "1"
        C = 1./float(Alpha) * Kfudge
        Time = float(data[1])
        RunTime[Run] = Time

        if Alpha != '1':
            if Run not in RunC:
                RunC[Run] = set()
            if C not in CRun:
                CRun[C] = set()
            RunC[Run].add(C)
            CRun[C].add(Run)


    for k, v in RunC.iteritems():
        if len(v) > 1:
            print "EXIT: More than one weight found for a single run!"
            exit(1)
    return RunC, CRun, RunTime


def FillOnOff(EvListPath, RunCArg):
    """ Returns OnOff[(C, EBin)][AreaType] where AreaType is either 'On' or 'Off' """
    OnOff = dict()
    for data in constants.get_event_line(EvListPath, ' '):
        Run = int(data[0])
        C = list(RunCArg[Run])[0]
        Alpha = data[9]
        if Alpha == "0":
            Alpha = "1"
        Energy = float(data[7])

        EBin = -1
        for iE in range(1, constants.NE):
            E1 = constants.Erange[iE]
            if Energy <= E1:
                EBin = iE - 1
                break

        if (C, EBin) not in OnOff:
            OnOff[(C, EBin)] = {'On': 0, 'Off': 0}

        if Alpha != '1':
            AreaType = 'Off'
        else: # Signal
            AreaType = 'On'
            # print "AreaType = On"

        OnOff[(C, EBin)][AreaType] = OnOff[(C, EBin)][AreaType] + 1
    return OnOff

def FillOnOff_EOnly(EvListPath):
    """ Returns OnOff[EBin][AreaType] where AreaType is either 'On' or 'Off' """
    OnOff = dict()
    for data in constants.get_event_line(EvListPath, ' '):
        Run = int(data[0])
        Alpha = data[9]
        if Alpha == "0":
            Alpha = "1"
        Energy = float(data[7])

        EBin = -1
        for iE in range(1, constants.NE):
            E1 = constants.Erange[iE]
            if Energy <= E1:
                EBin = iE - 1
                break

        if EBin not in OnOff:
            OnOff[EBin] = {'On': 0, 'Off': 0}

        if Alpha != '1':
            AreaType = 'Off'
        else: # Signal
            AreaType = 'On'
            # print "AreaType = On"

        OnOff[EBin][AreaType] = OnOff[EBin][AreaType] + 1
    return OnOff

def SaveOnOffPairs(OUT_DIR_ARG, Source, OnOffArg, AreasInfoArg, CRunArg, RunTimeArg):
    """ Save (On, Off) for each (C, EBin) in OnOffArg in form of 'obseration' objects. """
    for k, v in OnOffArg.iteritems():
        iE = k[1]
        if iE == -1:
            continue
        C         = k[0]
        Runs      = list(CRunArg[C])
        E0        = constants.Erange[iE]
        E1        = constants.Erange[iE + 1]
        Em        = (E0 + E1) / 2.
        TimesList = Times(C, CRunArg, RunTimeArg)
        # print "TimesList: ", TimesList
        TotalTime = sum(TimesList)
        AreaList  = Area(Em, Runs, AreasInfoArg)
        #print k, v, Runs, TimesList, AreaList, sum(map(lambda x, y: x * y, TimesList, AreaList))/sum(TimesList)
        MeanArea  = sum(map(lambda x, y: x * y, TimesList, AreaList))/TotalTime
        # print k, v, TotalTime, MeanArea
        print "k: ", k, "v: ", v, "Time: ", TotalTime, "Area: ", MeanArea
        if v['On'] + v['Off'] < 10:
            continue
        obs                  = bayes.Observation()
        obs.Ci               = C
        obs.Ti               = TotalTime
        obs.EnergyBin        = [E0, E1]
        obs.EnergyBinIndex   = iE
        obs.EnergyBinIndexes = [iE]
        obs.Npi              = v['On']
        obs.Nmi              = v['Off']
        obs.Ai               = [[MeanArea]]
        obs.P_Ai             = [[1.0]]
        obs.Ei               = [Em]
        obs.Eirange          = [(E0,E1)]
        obs.P_Ei             = [1.0]
        obs.source           = Source
        obs.Ji               = [constants.ALEX_JS_CONV_PSF["f"][Source][iE]]
        # obs.Ji               = [constants.ALEX_JS_CONV_PSF["f"][Source][30]]
        # obs.Ji               = [constants.J[Source]]
        obs.P_Ji             = [1.0]
        obs.save("%s/%s_%s_%s.pickle"\
            % (OUT_DIR_ARG, obs.EnergyBinIndex, obs.Ci, obs.source))

# Computations for Jim:
if False:
    OnOff_EOnly = FillOnOff_EOnly(PATHS["Segue"]["EvList"])
    print OnOff_EOnly
    AreasInfo = GetAEff(PATHS["Segue"]["Area"])
    RunC, CRun, RunTime = GetRunsCAndTime(PATHS["Segue"]["EvList"])
    AllSegRuns = [k for k, v in RunC.iteritems()]
    AllSegC    = map(lambda run: list(RunC[run])[0], AllSegRuns)
    print RunC

    TimesList = map(lambda x: RunTime[x], AllSegRuns)
    TotalTime = sum(TimesList)
    print "TotalTime: ", TotalTime
    MeanC  = sum(map(lambda x, y: x * y, TimesList, AllSegC))/TotalTime
    print "MeanC: ", MeanC
    for iE in range(constants.NE - 1):
        E0        = constants.Erange[iE]
        E1        = constants.Erange[iE + 1]
        Em        = (E0 + E1) / 2.
        AreaList  = Area(Em, AllSegRuns, AreasInfo)
        MeanArea  = sum(map(lambda x, y: x * y, TimesList, AreaList))/TotalTime
        print  iE, constants.Erange[iE], constants.Erange[iE + 1], MeanArea, OnOff_EOnly[iE]['On'], OnOff_EOnly[iE]['Off']
    exit(0)


for SourceInd in ["Segue", "Bootes1", "Draco", "UrsaMinor", "WilmanI"]:
    AreasInfo = GetAEff(PATHS[SourceInd]["Area"])
    RunC, CRun, RunTime = GetRunsCAndTime(PATHS[SourceInd]["EvList"])
    print "RunTime: ", RunTime
    OnOff = FillOnOff(PATHS[SourceInd]["EvList"], RunC)
    SaveOnOffPairs(OUT_DIR, SourceInd, OnOff, AreasInfo, CRun, RunTime)
