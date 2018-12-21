//-*-mode:c++; mode:font-lock;-*-
/**

* \class VAEventStats
* \ingroup Stage 6
* \brief Event statistics to aid lightcurve generation.
*
* Original Author: David Steele
*
* $Author: akuznetl $
* $Date: 2013/08/01 05:16:38 $
* $Revision: 1.18 $
* $Tag$
*
**/

#include <TFile.h>
#include "VAEventStats.h"

using namespace std;

VAEventStats::VAEventStats(TFile* f)
{
	//fSaveListsToASCII = false;
	
	// Set up event list files
	//fEventNumberFilenameOn = "eventListOn.txt";
	//fEventNumberFilenameOff = "eventListOff.txt";
	//fEventListOn.open( fEventNumberFilenameOn.c_str() );
	//fEventListOff.open( fEventNumberFilenameOff.c_str() );
	
	// Create root tree
	f->cd();
	pfEventStatsTree = new TTree("EventStatsTree", "Eventwise Statistics for Lightcurves");
	pfEventStatsTree->Branch("ArrayEventNum", &faArrayEventNum, "faArrayEventNum/i");
	pfEventStatsTree->Branch("RunNum", &faRunNum, "faRunNum/i");
	pfEventStatsTree->Branch("MJDInt", &faMJDInt, "faMJDInt/i");
	pfEventStatsTree->Branch("MJDDbl", &faMJDDbl, "faMJDDbl/D");
	pfEventStatsTree->Branch("DayNSDbl", &faDayNSDbl, "faDayNSDbl/D");
	pfEventStatsTree->Branch("OnEvent", &faOnEvent, "faOnEvent/O");
	pfEventStatsTree->Branch("OffEvent", &faOffEvent, "faOffEvent/O");
	pfEventStatsTree->Branch("Psi", &faPsi, "faPsi/D");
	pfEventStatsTree->Branch("Weight", &faEventWeight, "faEventWeight/D");
	pfEventStatsTree->Branch("Azimuth", &faAzimuthRad, "faAzimuthRad/F");
	pfEventStatsTree->Branch("Elevation", &faElevationRad, "faElevationRad/F");
	pfEventStatsTree->Branch("TrackingAzimuth", &faArrayTrackingAzimuth_Deg, "faArrayTrackingAzimuth_Deg/F");
	pfEventStatsTree->Branch("TrackingElevation", &faArrayTrackingElevation_Deg, "faArrayTrackingElevation_Deg/F");
	pfEventStatsTree->Branch("EnergyGeV", &faEnergy_GeV, "faEnergy_GeV/F");
	pfEventStatsTree->Branch("EnergyRMS", &faEnergyRMS_GeV, "faEnergyRMS_GeV/F");
	pfEventStatsTree->Branch("Noise", &faNoise, "faNoise/D");
	pfEventStatsTree->Branch("Offset", &faOffset, "faOffset/F");
	pfEventStatsTree->Branch("EffectiveArea", &faEffectiveArea, "faEffectiveArea/F");
	pfEventStatsTree->Branch("ThetaSq", &fThetaSq_Deg2, "fThetaSq_Deg2/F");
	pfEventStatsTree->Branch("PulsarPhase", &faPulsarPhase, "faPulsarPhase/D");
	pfEventStatsTree->Branch("MJDDblBarycentered", &faMJDDblBarycentered, "faMJDDblBarycentered/D");
	// BJZ Hack
	pfEventStatsTree->Branch("BDTScore", &faBDTScore, "faBDTScore/D");
}


VAEventStats::~VAEventStats()
{
}

void VAEventStats::fill(VAShowerData* pShower, bool fOnEvent, bool fOffEvent,  double fWeight, double noise, Float_t effectiveArea, Float_t theta2, Double_t phase, long double MJDBary, float offset, double fBDTScore)
{
	faOnEvent = fOnEvent;
	faOffEvent = fOffEvent;
	faEventWeight = (Double_t)fWeight;
	faPsi = sqrt(pShower->fDirectionXCamPlane_Deg * pShower->fDirectionXCamPlane_Deg
				 + pShower->fDirectionYCamPlane_Deg * pShower->fDirectionYCamPlane_Deg);
	faArrayEventNum = pShower->fArrayEventNum;
	faRunNum = pShower->fRunNum;
	faMJDInt = pShower->fTime.getMJDInt();
	faMJDDbl = pShower->fTime.getMJDDbl();
	faDayNSDbl = (Double_t)(1.*pShower->fTime.getDayNS());
	faEnergy_GeV = pShower->fEnergy_GeV;
	faEnergyRMS_GeV = pShower->fEnergyRMS_GeV;
	faElevationRad = pShower->fDirectionElevation_Rad;
	faAzimuthRad = pShower->fDirectionAzimuth_Rad;
	faArrayTrackingAzimuth_Deg = pShower->fArrayTrackingAzimuth_Deg;
	faArrayTrackingElevation_Deg = pShower->fArrayTrackingElevation_Deg;
	faNoise = noise;
	faEffectiveArea = effectiveArea;
	fThetaSq_Deg2 = theta2;
	faOffset = offset;
	
	faPulsarPhase = phase;
	faMJDDblBarycentered = MJDBary;

	faBDTScore = fBDTScore;
	
	pfEventStatsTree->Fill();
	
}


void VAEventStats::writeToFile(TFile* f)
{
	f->cd();
	pfEventStatsTree->Write();
}
