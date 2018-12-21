//-*-mode:c++; mode:font-lock;-*-
/*** \class vaStage6ssupp /data/lucifer1.1/veritas/lookups/vegas/v2_5_0/ea/ea_Oct2012_ua_ATM22_vegasv250rc5_7sam_050off_s700t2_std_MSW1p1_MSL1p3_MH7_ThetaSq0p01_LZA.root
 * \ingroup results extractor
 * \brief This is the main program of stage 6, which is designed
 *        for the visualization of the analysis results. Here is
 *        a short list of what it is currently doing:
 *        - applying the analysis cuts.
 *        - 1d single telescope analysis.
 *        - 2d single telescope analysis.
 *        - wobble analysis of the stereo array data.
 *        - ring background analysis of the stereo array data.
 *        it will include also the spectral evaluation procedure.
 *
 * Original author: Alexander Konopelko
 *
 * $Author: nahee $
 * $Date: 2014/09/10 20:30:37 $
 * $Revision: 1.79 $
 * $Tag$
 *
 **/
//
// Written by:
// A. Konopelko
// Physics Dept.
// Purdue Univ.
// West Lafayette, In. 479096
// akonopel@physics.purdue.edu
// 765-494-5394

// 28-November-2005
//
// c++ common headers.

#include <errno.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "VAStarCat.h"
#include "VAStage6Analysis.h"

using namespace std;

string VAStage6Analysis::sDefaultConfigDirectory = "config";
string VAStage6Analysis::sDefaultObsMode = "Wobble";
string VAStage6Analysis::sDefaultOutputFileName = "results_";
string VAStage6Analysis::sDefaultLookupFileName = "";
string VAStage6Analysis::sDefaultAcceptPlotLookup = "VAAcceptanceLookupPlot.root";
bool VAStage6Analysis::sDefaultDisplayStereoParamFlag = false;
bool VAStage6Analysis::sDefaultLightCurve = false;
double VAStage6Analysis::sDefaultRingSize = 0.13;
double VAStage6Analysis::sDefaultCBGRingUpper = -1.0;
double VAStage6Analysis::sDefaultCBGRingLower = -1.0;
bool VAStage6Analysis::sDefaultExcludeStars = true;
string VAStage6Analysis::sDefaultStarExclusionCatalog = "";
bool VAStage6Analysis::sDefaultExcludeSource = false;
double VAStage6Analysis::sDefaultStarMagLimit = 6.0;
bool VAStage6Analysis::sDefaultExcludeWobbleAnl = false;
bool VAStage6Analysis::sDefaultExcludeCBGAnl = true;
bool VAStage6Analysis::sDefaultExcludeRBMAnl = false;
bool VAStage6Analysis::sDefaultExcludeRBMFinal = false;
bool VAStage6Analysis::sDefaultExcludeRunStats = false;
bool VAStage6Analysis::sDefaultExcludeEventStats = false;
bool VAStage6Analysis::sDefaultSpectrum = false;
bool VAStage6Analysis::sDefaultUpperLimit = false;
//bool VAStage6Analysis::sDefaultEventOrBg = true;
string VAStage6Analysis::sDefaultBgEstimate = "wobble"; // Default spectrum analysis uses reflected ring
double VAStage6Analysis::sDefaultTestPositionRA = 0.;
double VAStage6Analysis::sDefaultTestPositionDEC = 0.;
double VAStage6Analysis::sDefaultStarExclusionRegionRadius = 0.0; // Default of 0 means use mag-dependent radii
double VAStage6Analysis::sDefaultSourceExclusionRegionRadius = 0.3; // [deg]
double VAStage6Analysis::sDefaultStarCatFilterRadius = 3.0;  //deg
string VAStage6Analysis::sDefaultUserExclusionRegionsFile = "";
int VAStage6Analysis::sDefaultDrawExclusionRegions = 0;
bool VAStage6Analysis::sDefaultPulsarAnlFlag = false;
bool VAStage6Analysis::sDefaultDoRelativeExposure = false;
int VAStage6Analysis::sDefaultNumRings = 0;
bool VAStage6Analysis::sDefaultReadFromStage4 = false;
bool VAStage6Analysis::sDefaultReadFromStage5Combined = false;
bool VAStage6Analysis::sDefaultBatchMode = false;
bool VAStage6Analysis::sDefaultSaveAcceptance = false;
bool VAStage6Analysis::sDefaultLoadAcceptance = false;
bool VAStage6Analysis::sDefaultApplyEnergyCorrection = false;
string VAStage6Analysis::sDefaultEnergyCorrectionFormula = "x";

VAStage6Analysis::VAStage6Analysis(int& argc, char** argv)
{

  VARegisterThisFunction vreg(__PRETTY_FUNCTION__);
  pfDebug =   VAGlobal::instance()->getDebugModule("VAStage6Analysis");
	
	
  fArgc = argc;
  fArgv = argv;
	
  fConfigDirectory = sDefaultConfigDirectory;
  fObsMode = sDefaultObsMode;
  fOutputFileName = sDefaultOutputFileName;
  fLookupFileName = sDefaultLookupFileName;
  fAcceptPlotLookup = sDefaultAcceptPlotLookup;
  fDisplayStereoParamFlag = sDefaultDisplayStereoParamFlag;
  fLightCurve = sDefaultLightCurve;
  fMCData = VAGlobal::instance()->isSimulationMode();
  fRingSize = sDefaultRingSize;
  fCBGRingUpper = sDefaultCBGRingUpper;
  fCBGRingLower = sDefaultCBGRingLower;
  fExcludeStars = sDefaultExcludeStars;
  fStarExclusionCatalog = sDefaultStarExclusionCatalog;
  fStarCatFilterRadius = sDefaultStarCatFilterRadius;
  fExcludeSource = sDefaultExcludeSource;
  fStarMagLimit = sDefaultStarMagLimit;
  fStarExclusionRegionRadius = sDefaultStarExclusionRegionRadius;
  fSourceExclusionRegionRadius = sDefaultSourceExclusionRegionRadius;
  fDrawExclusionRegions = sDefaultDrawExclusionRegions;
  fUserExclusionRegionsFile = sDefaultUserExclusionRegionsFile;
  fExcludeWobbleAnl = sDefaultExcludeWobbleAnl;
  fExcludeCBGAnl = sDefaultExcludeCBGAnl;
  fExcludeRBMAnl = sDefaultExcludeRBMAnl;
  fExcludeRBMFinal = sDefaultExcludeRBMFinal;
  fExcludeRunStats = sDefaultExcludeRunStats;
  fExcludeEventStats = sDefaultExcludeEventStats;
  fSpectrum = sDefaultSpectrum;
  fUpperLimit = sDefaultUpperLimit;
  //fEventOrBg = sDefaultEventOrBg;
  fBgEstimate = sDefaultBgEstimate;
  fTestPositionRA = sDefaultTestPositionRA;
  fTestPositionDEC = sDefaultTestPositionDEC;
  fPulsarAnlFlag = sDefaultPulsarAnlFlag;
  fDoRelativeExposure = sDefaultDoRelativeExposure;
  fNumRings = sDefaultNumRings;
  fReadFromStage4 = sDefaultReadFromStage4;
  fReadFromStage5Combined = sDefaultReadFromStage5Combined;
  fBatchMode = sDefaultBatchMode;
  fLoadAcceptance = sDefaultLoadAcceptance;
  fSaveAcceptance = sDefaultSaveAcceptance;
  fApplyEnergyCorrection = sDefaultApplyEnergyCorrection;
  fEnergyCorrectionFormula = sDefaultEnergyCorrectionFormula;
	
	
	
  fSourceName = "Source";
	
	
  fCutMask = 0;
  fTestPosition = false;
  fTestRA = 0.;
  fTestDec = 0.;
  fSourceRA = 0;
  fSourceDec = 0;
	
  fExposureOn = 1.e-10;
  fExposureOff = 1.e-10;
  fLineCount = 0;
  fErrstr.str("");
	
  fRun = "On";
	
  pfStereoCutter = new VAStereoCuts;
  pfHillasCutter = NULL;
  pfTelMaskCutter = VATelMaskCutterFactory::instance()->
    getTelMaskCutter(pfHillasCutter);
  pfTelMaskCutter->reset();
	
  pfSizeMaskCutter = VASizeMaskCutterFactory::instance()->
    getSizeMaskCutter();
  pfSizeMaskCutter->reset();
	
  pfSizeMaskCutter->printAllowedMasks();
	
  pfLightCurveAnl = new VALightCurveAnl();
  pfShowerParamDisplay = new VAShowerParamDisplay();
  pfPulsar = new VAPulsarAnl();
	
  VAEffectiveAreaManager* pfEAManagerTemp  = new VAEffectiveAreaManager();
  pfMigrationMatrix = (TH2F*)pfEAManagerTemp->getMigrationMatrix()->Clone();
	
  delete pfEAManagerTemp;
  pfEAParameterData = new VAEAParameterData(); // Needed by the EA manager code
	
  fRunNumber = 0;
  fNgroup = 0;
  fGroupID.clear();
}

void VAStage6Analysis::configure(ConfigInfo& config_file, VSOptions& command_line)
{
  // set default config directory.
  doVAConfiguration(config_file, command_line,
		    "S6A_ConfigDir", sDefaultConfigDirectory,
		    "Stage 6 Main",
		    "Here you should specify the configuration directory. "
		    "If you want to save the configuration file by default "
		    "you have to specify configuration directory first. All "
		    "output files of stage 6 will be saved in this directory. "
		    "It is strongly recommended that you create the configuration "
		    "directory before running stage6 code.", true);
					  
  // set the name for the output file.
  doVAConfiguration(config_file, command_line,
		    "S6A_OutputFileName", sDefaultOutputFileName, "Stage 6 Main",
		    "All histograms and plots will be saved into the output root "
		    "file with this name. In addition, the plots will be printed "
		    "into the postscript file with the same name too. ", true);
					  
  // set the observational mode.
  doVAConfiguration(config_file, command_line,
		    "S6A_ObsMode", sDefaultObsMode, "Stage 6 Main",
		    "There are two observational mode used so far with VERITAS: "
		    "On/Off and Wobble. If you chose the On/Off mode the run    "
		    "list should be composed in a sequence of On and Off runs, "
		    "respectively.  Select Wobble for all observations other "
		    "than On/Off observations.", false);
					  
					  
  // set the flag for displaying stereo parameters.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_StereoParamDisplay", sDefaultDisplayStereoParamFlag,
		    "Stage 6 Main",
		    "If you want to display a number of generic parameters "
		    "used for the stereoscopic analysis e.g. mean scaled Width "
		    "distribution of the impact distances etc, please set this "
		    "key to true or 1.", false);
					  
  // angular size of the rings in the Crescent mode analysis.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_CBGRingUpper", sDefaultCBGRingUpper, "Stage 6 Main",
		    "This defines the angular size of the upper ring used in the crescent background model. ", -1.0);
					  
  doVAConfiguration(config_file, command_line,
		    "S6A_CBGRingLower", sDefaultCBGRingLower, "Stage 6 Main",
		    "This defines the angular size of the lower ring used in the crescent background model. ", -1.0);
					  
  // angular size of the ring in the Wobble mode analysis.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_RingSize", sDefaultRingSize, "Stage 6 Main",
		    "This defines the angular size of the regions used in the reflected region model. "
		    "Note that this value is the square root of the quantity often referred to as theta squared", true);
					  
  doVAConfiguration(config_file, command_line,
		    "S6A_NumRings", sDefaultNumRings, "Stage 6 Main",
		    "This defines the number of integration rings used in the reflected region model.", true);
					  
  // flag for use of exclusion regions for stars.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_ExcludeStars", sDefaultExcludeStars,
		    "Stage 6 Main",
		    "This defines whether star exclusion regions "
		    "are used.  Star exclusions are on (true) by default.", true);
					  
  // set the name for the collection areas file.
  doVAConfiguration(config_file, command_line,
		    "S6A_LookupFileName", sDefaultLookupFileName,
		    "Stage 6 Main",
		    "For evaluation of the energy spectrum an auxilary file "
		    "which contains the effective collection areas, is needed."
		    "This parameter defines the name of this file. ", true);
					  
  // flag for use of exclusion regions for stars.
  doVAConfiguration(config_file, command_line,
		    "S6A_StarExclusionCatalog", sDefaultStarExclusionCatalog,
		    "Stage 6 Main",
		    "By default, the Tycho catalog is used, and read from vegas/resultsExtractor/utilities. "
		    "If you need to read a star catalog from a different file or location, use this flag to override "
		    "the default location (should only need to do so under special circumstances, like using a distributed computing network). "
		    "Note that this is for stars; for user-defined exclusions, use S6A_UserDefinedExclusionList.", true);
					  
  // flag to set the radius of stars included for exclusion
	
  doVAConfiguration(config_file, command_line,
		    "S6A_StarCatFilterRadius", sDefaultStarCatFilterRadius,
		    "Stage 6 Main",
		    "Use this to choose the radius of stars to consider for exclusion regions. "
		    "Set to <= 0 for NO RADIUS FILTERING.  Default is 3.0 degrees.", true);
					  
					  
  // different options for drawing exclusion regions and stars on the skymaps
  doVAConfiguration(config_file, command_line,
		    "S6A_DrawExclusionRegions", sDefaultDrawExclusionRegions,
		    "Stage 6 Main",
		    "Use this to draw the exclusion regions on the skymaps. "
		    "Options are: 0=none, 1=regions alone, 2=stars alone, 3=both. Default is 0", true);
					  
  // flag for use of exclusion regions for putative source.
  doVAConfiguration(config_file, command_line,
		    "S6A_ExcludeSource", sDefaultExcludeSource,
		    "Stage 6 Main",
		    "This defines whether an exclusion region for the "
		    "source is used - source taken from run header of "
		    "first run, or assumed to be the TestPosition coordinates the the TestPositionRA and TestPositionDec flags are used.  Set true to use "
		    "exclusion region.  Note: runs taken in Survey mode define the 'source' as the tracking direction so with them, "
		    "this flag should only be used if you are also using the TestPosition flags.  See also S6A_UserDefinedExclusionList.", true);
					  
  // Angular radius of the star exclusion region.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_StarExclusionRadius", sDefaultStarExclusionRegionRadius,
		    "Stage 6 Main",
		    "This is the radius of the STAR exclusion region in deg."
		    "If you set this, then ONE value is used for all stars, independent of magnitude. "
		    "The default is a 'smart' algorithm based on magnitude with 0.3 for most stars.", true);
					  
  // Angular radius of the exclusion region.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SourceExclusionRadius", sDefaultSourceExclusionRegionRadius,
		    "Stage 6 Main",
		    "This is the radius of the SOURCE exclusion region in deg. The default 0.3 seems reasonable for moderately bright point sources; "
		    "0.4 better for Crab-strength point sources.", true);
					  
  // Name of the file containing user-defined exclusion region infos.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_UserDefinedExclusionList", sDefaultUserExclusionRegionsFile , "Stage 6 Main",
		    "Here you define the filename of a list of sources that should be excluded from background estimation. "
		    "Each row in the file defines one exclusion and requires four fields: RA and Dec in decimal degrees, "
		    "radius of exclusion in degrees, and a name for the exclusion. "
		    "A list of stars to exclude can be generated by using the S6A_ExcludeStars=true flag.", true);
					  
  // Magnitude limit for excluding stars (B-band).
	
  doVAConfiguration(config_file, command_line,
		    "S6A_StarExclusionBMagLimit", sDefaultStarMagLimit,
		    "Stage 6 Main",
		    "Sets the magnitude cutoff (B-band) used to generate a list of stars to exclude from a field."
		    "Note that the optimal value has not yet been carefully studied: the default is set to include stars of mag < 6.0 based "
		    "on prior experience, but this may yield too many exclusions in galactic fields, "
		    "for which a limit of 5.0 or 5.5 may be more appropriate.  Be wary!", true);
					  
  // flag to suppress the wobble mode computations.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressWobble", sDefaultExcludeWobbleAnl,
		    "Stage 6 Main",
		    "For some sky survey runs it makes sense to supress "
		    "the final Wobble Analysis computations.", false);
					  
  // flag to suppress the crescent background model computations.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressCBG", sDefaultExcludeCBGAnl,
		    "Stage 6 Main", "Can now choose to run Crescent analysis instead of Wobble.", true);
					  
  // set the name for the acceptance plot lookup file.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_AcceptanceLookup", sDefaultAcceptPlotLookup, "Stage 6 Main",
		    "The radial acceptance plot can be saved or replaced."
		    "This is the name of the data file to which the "
		    "acceptance is written or from which it is read, if "
		    "either the saveAcceptance or loadAcceptance flags "
		    "are set.", false);
					  
  doVAConfiguration(config_file, command_line,
		    "S6A_LoadAcceptance", sDefaultLoadAcceptance, "Stage 6 Main",
		    "Load acceptance curves. Goes with kAcceptanceLookup=true", false);
					  
  doVAConfiguration(config_file, command_line,
		    "S6A_SaveAcceptance", sDefaultSaveAcceptance, "Stage 6 Main",
		    "Save acceptance curves. Goes with kAcceptanceLookup=true", false);
					  
					  
  // flag to suppress the Ring Background Model computations.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressRBM", sDefaultExcludeRBMAnl,
		    "Stage 6 Main",
		    "Sometimes the RBM computations can take a lot of time."
		    "If you do not need these results, RBM can be turned "
		    "off by setting this to true.", true);
					  
  // RA of the 'test position'.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_TestPositionRA", sDefaultTestPositionRA,
		    "Stage 6 Main",
		    "Here you define the RA coordinate of the center "
		    "of the sky maps defined in RBM, as well as the RA "
		    "of the ON region for the Wobble Analysis. Be aware "
		    "that "
		    "the spectrum will be computed for this source then. "
			"Units are in degrees.", true);
					  
  // DEC of the 'test position'.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_TestPositionDEC", sDefaultTestPositionDEC,
		    "Stage 6 Main",
		    "Here you define the DEC coordinate of the center "
		    "of the sky maps defined in RBM, as well as the Dec of "
		    " the ON region for the Wobble Analysis. Be aware that "
		    "the spectrum will be computed for this source then. "
			"Units are in degrees.", true);
					  
  // flag to suppress just the final stage of RBM computation.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressRBMFinalStage", sDefaultExcludeRBMFinal,
		    "Stage 6 Main",
		    "This flag suppresses only the final stage (map "
		    "generation - the long part - of the Ring Background "
		    "Model.  Use this if you want to suppress the RBM "
		    "map generation but plan to use the relative "
		    "exposure corrections (kDoRelativeExposure=1)", false);
					  
  // flag to suppress the VARunStats computations.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressRunStats", sDefaultExcludeRunStats,
		    "Stage 6 Main",
		    "By default RunStats (runwise descriptive statistics"
		    " are generated (0), but they can be suppressed by"
		    " setting this flag to true (1).", false);
					  
  // flag to suppress the VAEventStats writing
	
  doVAConfiguration(config_file, command_line,
		    "S6A_SuppressEventStats", sDefaultExcludeEventStats,
		    "Stage 6 Main",
		    "By default EventStats (limited eventwise statistics"
		    " are not saved (1), but they can be saved in a root file by"
		    " setting this flag to false (0).", false);
					  
					  
  // if user selected the spectrum more, the spectrum will be generated.
  doVAConfiguration(config_file, command_line,
		    "S6A_Spectrum", sDefaultSpectrum,
		    "Stage 6 Main",
		    "This defines whether the spectrum will be "
		    "calculated or not. Select -true- if you want "
		    "the spectrum to be calculated.", true);
					  
  // if user wants an upper limit.
  doVAConfiguration(config_file, command_line,
		    "S6A_UpperLimit", sDefaultUpperLimit,
		    "Stage 6 Main",
		    "This defines whether an upper limit will be "
		    "calculated or not. Select -true- if you want "
		    "the upper limit to be calculated. This is independent on the spectrum, "
		    "that means that you may calculate both a spectrum and an upper limit over the whole "
		    "data set. This is NOT the option to set the calculation of an upper limit "
		    "beyond the last significant spectrum point, use SP_UL instead.", false);
					  
					  
					  
  /*	// this define which class to use for events extraction while doing spectrum.
	
	doVAConfiguration(config_file, command_line,
	"S6A_BgEstimate", sDefaultEventOrBg,
	"Stage 6 Main",
	"This defines whether the spectrum will use "
	"the events selected by the RBM (0) or Wobble "
	"Analysis (1, default).", false);
  */
  // this define which class to use for events extraction while doing spectrum.
	
  doVAConfiguration(config_file, command_line, "S6A_BgEstimate", sDefaultBgEstimate,
		    "Stage 6 Main", "This defines whether the spectrum "
		    "will use the events selected by the RBM (ring), Wobble "
		    "(wobble), or CBG (crescent) analysis (wobble, default).", "wobble");
					  
  // set the flag for the light curve
	
  doVAConfiguration(config_file, command_line,
		    "S6A_LightCurve", sDefaultLightCurve, "Stage 6 Main",
		    "To generate the Light Curve please enable this key. "
		    "Please remember that the right time sequence of the input "
		    "data runs is important here! This light-curve code is "
		    "currently broken; do not use.", false);
					  
					  
					  
  // set the flag for pulsar analysis.
	
  doVAConfiguration(config_file, command_line,
		    "S6A_PulsarAnl", sDefaultPulsarAnlFlag, "Stage 6 Main",
		    "This flag allows you to do the Pulsar Analysis. ", false);
					  
					  
					  
					  
  // set the flag for the relative exposure calculation.
  doVAConfiguration(config_file, command_line,
		    "S6A_DoRelativeExposure", sDefaultDoRelativeExposure,
		    "Stage 6 Main",
		    "This flag enables calculation of relative exposure "
		    "between runs, which takes care of weighting runs "
		    "properly in wobble, light curve, spectrum when alpha "
		    "varies between runs.  Set to false by default while "
		    "under development.", false);
					  
  // set the flag for the relative exposure calculation.
  doVAConfiguration(config_file, command_line,
		    "S6A_ReadFromStage4", sDefaultReadFromStage4,
		    "Stage 6 Main",
		    "Set to true if you want to read shower data directly from stage 4, without running stage 5. "
		    "It is recommended that you do run stage 5, as the analysis is faster in the long run", true);
					  
  // set the flag for the relative exposure calculation.
  doVAConfiguration(config_file, command_line,
		    "S6A_ReadFromStage5Combined", sDefaultReadFromStage5Combined,
		    "Stage 6 Main",
		    "Set to true if you want to read shower data and other information from the Stage5 combined tree "
		    "This may become the recommended default", true);
					  
  // set the flag for the relative exposure calculation.
  doVAConfiguration(config_file, command_line,
		    "S6A_Batch", sDefaultBatchMode,
		    "Stage 6 Main",
		    "Set this to true if you don't want any plots and want stage 6 to "
		    "return to the command line upon completion. Useful for large "
		    "batch jobs for eg optimisation", true);
					  
  // set flag to enable correction for the energy scale
  doVAConfiguration(config_file, command_line,
		    "S6A_ApplyEnergyCorrection", sDefaultApplyEnergyCorrection,
		    "Stage 6 Main",
		    "Set this flag if you want to apply an energy scale correction. "
		    "The formula to modify the energy scale is defined in "
		    "S6A_EnergyCorrectionFormula option. Default false. "
		    "ATTENZIONE!!!!! "
		    "If you decide to modify the energy scale, your spectrum "
		    "and upper limit will be modified accordingly, "
		    "therefore will be written in the \"IndecentResults\" directory, "
		    "in order to warn you that the result has been altered.", true);
					  
  // set the energy scale modification formula
  doVAConfiguration(config_file, command_line,
		    "S6A_EnergyCorrectionFormula", sDefaultEnergyCorrectionFormula,
		    "Stage 6 Main",
		    "Set a string with the formula used to modify the energy scale. "
		    "The format for the formula is the same as for a TFormula object. "
		    "This option is ignored if S6A_ApplyEnergyCorrection is set to false. "
		    "Default \"x\" (no modification of the energy scale).", true);
					  
}


bool VAStage6Analysis::setup(ConfigInfo* config, ConfigInfo* cuts)
{

  if(fDoRelativeExposure && fExcludeRBMAnl)
    {
      fExcludeRBMAnl = false;
      fExcludeRBMFinal = true;
      fErrstr.str("");
      fErrstr << "You chose to suppress the Ring Background Model "
	      << "completely at the same time as using the relative "
	      << "exposure calculation, which relies on the runwise "
	      << "acceptance curves generated by the RBM.  Let's make a "
	      << "deal: I'll skip the final, slow stage of the RBM but "
	      << "generate the acceptance curves for the relative exposure "
	      << "calculation.  In the future, you can suppress just the "
	      << "final stage of the RBM analysis with the flag "
	      << " -kSuppressRBMFinalStage=1.";
      VAGlobal::instance()->getErrorReport()
	->addError(fErrstr.str().c_str());
    }
	
  if(fSpectrum)
    {
      if(fBgEstimate == "crescent")
	{
	  if(fExcludeCBGAnl)
	    {
	      fExcludeCBGAnl = false;
	      fErrstr.str("");
	      fErrstr << "You chose to run the spectrum analysis using the Crescent "
		      << "analysis for background estimation (-S6A_BgEstimate=crescent) "
		      << "but disabled the crescent analysis (-S6A_SuppressCBG=1). "
		      << "I will override that and turn on the crescent analysis so "
		      << "the spectral analysis will work. Did you mean to use the "
		      << "Wobble analysis for background estimate (-S6A_BgEstimate=wobble) "
		      << "or perhaps the RBM (-S6A_BgEstimate=ring)?";
	      VAGlobal::instance()->getErrorReport()
		->addError(fErrstr.str().c_str());
	    }
	}
      else if(fBgEstimate == "wobble")
	{
	  if(fExcludeWobbleAnl)
	    {
	      fExcludeWobbleAnl = false;
	      fErrstr.str("");
	      fErrstr << "You chose to run the spectrum analysis using the Wobble "
		      << "analysis for background estimation (-S6A_BgEstimate=0) "
		      << "but disabled the wobble analysis (-S6A_SuppressWobble=1). "
		      << "I will override that and turn on the wobble analysis so "
		      << "the spectral analysis will work. Did you mean to use the "
		      << "CBG method for the background estimate (-S6A_BgEstimate=crescent) "
		      << "or perhaps the RBM analysis (-S6A_BgEstimate=ring)?";
	      VAGlobal::instance()->getErrorReport()
		->addError(fErrstr.str().c_str());
	    }
	}
      else if(fBgEstimate == "ring")
	{
	  if(fExcludeRBMAnl)
	    {
	      fExcludeRBMAnl = false;
	      fErrstr.str("");
	      fErrstr << "You chose to run the spectrum analysis using the RBM "
		      << "analysis for background estimation (-S6A_BgEstimate=ring) "
		      << "but disabled the RBM analysis (-S6A_SuppressRBM=1). "
		      << "I will override that and turn on the RBM analysis so "
		      << "the spectral analysis will work, but I will turn off the "
		      << "final stage of the RBM analysis "
		      << "(-S6A_SuppressRBMFinalStage=1) to skip the map "
		      << "generation. Did you mean to use the Wobble analysis "
		      << "(-S6A_BgEstimate=wobble) or CBG analysis (-S6A_BgEstimate=crescent) "
		      << "for background estimate?";
	      VAGlobal::instance()->getErrorReport()
		->addError(fErrstr.str().c_str());
	    }
	}
      else
	{
	  ostringstream errstr("");
	  errstr << "+++ S6A WARNING: You set the backgound mode improperly (-S6A_BgEstimate=?). Please set to wobble, crescent, or ring" << std::endl;
	  VAGlobal::instance()->getErrorReport()->addError(errstr.str().c_str());
	  std::cout << errstr.str().c_str();
	}
    }
	
  if((!fExcludeRBMAnl && !fExcludeWobbleAnl) || (!fExcludeCBGAnl && !fExcludeWobbleAnl) || (!fExcludeRBMAnl && !fExcludeCBGAnl && !fExcludeWobbleAnl) || (!fExcludeRBMAnl && !fExcludeCBGAnl))
    {
      ostringstream errstr("");
      errstr << "+++ S6A: WARNING: You set to do multiple analyses, in this case the UL will be calculated using the method specified by -S6A_BgEstimate (default is wobble). " << std::endl;
		
      errstr << "                  If you want to calculate the UL from another method, please set -S6A_BgEstimate=(ring, wobble, crescent)." << std::endl;
#ifndef __CINT__
      //VAGlobal::instance()->getErrorReport()->addError(errstr.str().c_str());
#endif
      std::cout << errstr.str().c_str() << endl;
		
    }
	
  // check out the config directory.
	
  DIR* pdir;
  struct dirent* pent;
  pdir = opendir(fConfigDirectory.c_str());
  if(!pdir)
    {
      fErrstr.str("");
      if(fConfigDirectory == "config")
	{
	  fErrstr << "+++ Default config directory : '" << fConfigDirectory
		  << "' does not exist!    " << endl;
	}
      else
	{
	  fErrstr << "+++ You specified config directory : '" << fConfigDirectory
		  << "' which does not exist!    " << endl;
	}
      fErrstr << "+++ Please create the directory first" << endl;
      fErrstr << "+++ mkdir " << fConfigDirectory << endl;
      fErrstr << "+++ in that place from where you run the code!" << endl;
      VAGlobal::instance()->getErrorReport()
	->addError(fErrstr.str().c_str());
      VAGlobal::instance()->getErrorReport()
	->printErrorReport();
      exit(EXIT_SUCCESS);
    }
  else
    {
      cout << "+++ Config directory : " << fConfigDirectory << endl;
    }
	
  errno = 0;
  ostringstream os("");
	
  while((pent = readdir(pdir)))
    {
      os << "+++ file : " << pent->d_name << endl;
    }
  pfDebug->DebugMessage(E_INFO, os.str().c_str());
  if(errno)
    {
      VAGlobal::instance()->getErrorReport()
	->addError("+++ Reading config dir failured, terminating!");
      VAGlobal::instance()->getErrorReport()
	->printErrorReport();
      exit(1);
    }
  closedir(pdir);
	
	
  pfStereoCutter->printCuts();
	
	
  // check the test position.
	
  double eps = FLT_EPSILON;
	
  if((fabs(fTestPositionRA) > eps) &&
     (fabs(fTestPositionDEC) > eps))
    {
      fTestPosition = true;
      fTestRA  = fTestPositionRA * TMath::Pi() / 180.;
      fTestDec = fTestPositionDEC * TMath::Pi() / 180.;
    }
  if(fTestPosition)
    {
      cout << "+++ All Sky Maps will be centered at the Test Position: " << endl;
      cout << "+++ RA : " << fTestRA << " DEC: " << fTestDec << endl;
    }
	
  // now get the run list file.
	
  fArgv++;
  fArgc--;
  if(fArgc == 0)
    {
      fErrstr.str("");
      fErrstr << "+++ No Input Run file was specified!" << endl;
      fErrstr << "+++ To do analysis you need to specify the run list"
	      << endl;
      fErrstr << "+++ Example: bin/vaStage6 [options] mkn421_stereo_runlist "
	      << endl;
      VAGlobal::instance()->getErrorReport()
	->addError(fErrstr.str().c_str());
      VAGlobal::instance()->getErrorReport()
	->printErrorReport();
      exit(EXIT_FAILURE);
    }
	
  // import the name of the input file with the run list.
	
  string FileListName = *fArgv;
  ifstream fin;
  string tmps;
  uint32_t pGroupID = 1;
  uint32_t NumOfGroup = 0;
  fin.open(FileListName.c_str(), ifstream::in); //open runlist file
  bool runlistFlag = true;
  uint32_t Nrunfile = 0;
	
  //nahee:let's make it a bit complicated here...
  if(!fin.fail())
    {
      bool inComment = false;
      char lastc = '\n';
      for(;;)
        {
          char c = fin.get();
          if(!fin.good())
            {
	      break;	// EOF
	    }

          if(c == '#')
            {
              inComment = true;
            }
          else if(c == '\n')
            {
              inComment = false;
            }
          if(!inComment && (lastc != '\n' || c > ' ') )
            {
              fRunlists << c;
              lastc = c;
            }
        }

      fin.close();
      while(getline(fRunlists, tmps)) //read runlist line by line
	{
	  if((tmps.substr(0, 1) == "[") && (tmps.substr(0, 2) != "[/") && tmps.find("RUNLIST") < tmps.size())
	    {
	      sscanf(tmps.c_str(), "[%*s %*s %d]", &pGroupID); //read ID #
	      for(unsigned int n = 0; n < fGroupID.size(); n++)
		{
		  if(pGroupID == fGroupID.at(n))
		    {
		      fErrstr.str("");
		      fErrstr << "+++ Apparently you used identical Group ID " << pGroupID
			      << " for multiple times." << endl;
		      fErrstr << "+++ Please check again." << endl;
		      VAGlobal::instance()->getErrorReport()
			->addError(fErrstr.str().c_str());
		      VAGlobal::instance()->getErrorReport()
			->printErrorReport();
		      exit(EXIT_SUCCESS);
		    }
		}
	      fGroupID.push_back(pGroupID);
	      NumOfGroup ++;
	      runlistFlag = true;
	    }
	  else if(tmps.substr(0, 2) == "[/" && tmps.find("RUNLIST") < tmps.size())
	    {
	      runlistFlag = false;
	    }
	  else if((tmps.substr(0, 1) == "[") && (tmps.substr(0, 2) != "[/") && tmps.find("RUNLIST") >= tmps.size() && runlistFlag)
	    {
	      runlistFlag = false;
	    }
	  else if(runlistFlag)
	    {
	      Nrunfile ++;
	    }
			
	  fLineCount++;
	  if(fLineCount == 1 && Nrunfile == 1)
	    {
	      //use default Group ID of 0 for first files
	      fGroupID.push_back(0);
	      NumOfGroup ++;
	    }
	} //end while loop (over runlist file)
    } //end if (!fin.fail())
  else
    {
      //fin fails
      fErrstr.str("");
      fErrstr << "+++ Apparently you specified a run list file"
	      << endl;
      fErrstr << "+++ which does not exist!"
	      << endl;
      fErrstr << "+++ Please run again with the correct file."
	      << endl;
      VAGlobal::instance()->getErrorReport()
	->addError(fErrstr.str().c_str());
      VAGlobal::instance()->getErrorReport()
	->printErrorReport();
      exit(EXIT_SUCCESS);
    } // end else
  if(NumOfGroup == 1 && fGroupID.size() == 1)
    {
      cout << "+++ Single group running! " << endl;
      fNgroup = 1;
    }
  else if(NumOfGroup != fGroupID.size())
    {
      cout << "+++ Check the group IDs. Found " << NumOfGroup << " & " << fGroupID.size() << endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      cout << "+++ Found " << NumOfGroup << " groups! " << endl;
      fNgroup = NumOfGroup;
    }
	
  // summary here
  cout << "+++ Number of lines in the input file '"
       << FileListName.c_str() << "' : "
       << fLineCount << endl;
  cout << "+++ Number of data files in the input file : "
       << Nrunfile << endl;
  string::size_type idx = FileListName.find(".root");
  if(idx != string::npos)
    {
      //find returns npos if .root extension is not found
      VAGlobal::instance()->getErrorReport()
	->addError("+++ Please submit correct name of run list file!!!");
      VAGlobal::instance()->getErrorReport()
	->printErrorReport();
      exit(EXIT_FAILURE);
    }
    
  // temporary coordinates for simulated events are fixed!!!
  if(fMCData)
    {
      //CORSIKA
      fSourceRA =  6.06271;
      fSourceDec = 0.55768;
      //KASCADE
      fSourceRA =  345.8 / 180.*3.14159265;
      fSourceDec = 32.6 / 180.*3.14159265;
    }
	
  // Initialize the output ROOT file.
  idx = fOutputFileName.find('.');
  if(idx != string::npos)
    {
      fOutputFileName = fOutputFileName.substr(0, idx);
    }
  string outfile = fConfigDirectory + '/' + fOutputFileName + "s6.root";
  cout << "+++ Results will be written to ./"
       << outfile << endl;
  pfRootIO = new VARootIO(outfile);
  if(!pfRootIO)
    {
      throw VAFatalException("Unable to create the VARootIO object!!", __FILE__, __LINE__);
    }
  pfRootIO->initTheRootFile(false);
  if(!pfRootIO->testFileOpenForIO())
    {
      throw VAFatalException("Unable to create a writeable VARootIO object!!", __FILE__, __LINE__);
    }
	
  pfOutputFile = NULL;
  pfOutputFile = pfRootIO->getTFilePtr();//new TFile(outfile.c_str(),"recreate");
  if(!pfOutputFile)
    {
      throw VAFatalException("Unable to get a pointer to the TFile!!", __FILE__, __LINE__);
    }
	
  ostringstream configfile_stream;
  config->dumpConfigFile(configfile_stream);
  ostringstream cutsfile_stream;
  cuts->dumpConfigFile(cutsfile_stream);
  pfRootIO->setSomeConfigFileText(6, VAGlobal::instance()->getVersion(), configfile_stream.str());
  pfRootIO->writeConfigFile();
  pfRootIO->setSomeCutsFileText(6, cutsfile_stream.str());
  pfRootIO->writeCutsFile();
	
	
  // Initialize the runwise statistics object
  pfDebug->DebugMessage(E_FULLDEBUG, "Calling VARunStats constructor");
	
  pfRunStats = new VARunStats(!fReadFromStage4, pfOutputFile);
  pfEventStats = NULL;
	
  if(!fExcludeEventStats)
    {
      // Initialize the eventwise statisitics object
      pfEventStats = new VAEventStats(pfOutputFile);
    }
	
  if(!fExcludeCBGAnl)
    {
      double theta2 = fRingSize * fRingSize;
      pfCrescentAnl = new VACrescentBackgroundAnl(pfRunStats, pfEventStats, fSourceExclusionRegionRadius, theta2);
    }
	
  pfHistogram2DAnl = new VAFill2DHistogramAnl(fExcludeSource, fSourceExclusionRegionRadius, pfOutputFile);
  pfWobbleAnl = new VAWobbleAnl(pfRunStats);
	
  pfSpectrum = new VASpectrumAnl(true, "Stage6Spectrum", "Spectrum", fMCData);
  if(fApplyEnergyCorrection)
    {
      pfSpectrum->SetName("VAIndecentSpectrum");
      //pfSpectrum->SetSpectrumName("IndecentSpectrum");
      //string literal  within " " are type const char*
      //gcc after v4.4 gives a warning if you try to do the above
      // so instead type cast it
      pfSpectrum->SetSpectrumName((char*)"IndecentSpectrum");
    }
  pfLiveTime = new VALiveTime(true);
  pfSpectrum->SetLiveTime(pfLiveTime);
  pfUpperLimit = new VAUpperLimit();
  pfUpperLimit->SetLiveTime(pfLiveTime);
	
  pfExposureCalc = new VAExposureCalculator(pfRunStats);
  if(fDoRelativeExposure &&
     !pfExposureCalc->getLibraryName().size())
    {
      fDoRelativeExposure = false;
      fErrstr.str("");
      fErrstr << "\033[33;1m Acceptance library file not specified!  Exposure "
	      << "calculation will ignore acceptance effects.\033[0m";
      VAGlobal::instance()->getErrorReport()
	->addError(fErrstr.str().c_str());
    }
  // read the run file from the list.
  cout << "+++ Read run name from a file list: "
       << FileListName << endl;
		 
  string RunFileName;
  // open just 1st file in the runlist.
  // load all avilable trees in the input file.
  fRunlists.clear();
  fRunlists.seekg(0, fRunlists.beg);
  fRunlists >> RunFileName;
	
  if(!fMCData) //is not monte carlo simulation
    {
      VARootIO* pRootIOtemp = new VARootIO(RunFileName, true); // true for read only
      //VAGlobal::instance()->addRootFile(&pRootIOtemp);
      pRootIOtemp->loadTheRootFile();
      if(!pRootIOtemp->IsOpen())
	{
	  throw VAFatalException("There was a problem loading file " + RunFileName, __FILE__, __LINE__);
	}
		
      VARunHeader* pRunHeader = pRootIOtemp->loadTheRunHeader();
      if(!pRunHeader)
	{
	  throw VAFatalException("There was a problem loading the run header from file " + RunFileName, __FILE__, __LINE__);
	}
      /*
      // We load the effective areas only if the user wants to do a spectral analysis
      if ((fSpectrum || fUpperLimit) && fNgroup == 1 ){
      VARootIO *pfEARootIO = new VARootIO(fLookupFileName, true);
      VAGlobal::instance()->addRootFile(&pfEARootIO);
      pfEARootIO->loadTheRootFile();
      if (!pfEARootIO->IsOpen())
      throw VAFatalException("Couldn't load the effective area tables from "+fLookupFileName, __FILE__, __LINE__);
      pfEAManager->loadEffectiveAreas(pfEARootIO);
      pfEAManager->getConfig().getValue("WindowForNoise", fWindowSizeForNoise);
      cout<<"Stage 6 Analysis picked up a window size for noise of "<<fWindowSizeForNoise<<endl;
      }
      */
      fSourceRA  = pRunHeader->fSourceInfo.fRA;
      fSourceDec = pRunHeader->fSourceInfo.fDec;
		
      delete pRunHeader;
      pRunHeader = NULL;
      pRootIOtemp->closeTheRootFile();
      delete pRootIOtemp;
      pRootIOtemp = NULL;
    }
	
	
  // the central point of the sky maps can be chosen by the user,
  // otherwise the central point coincides with the source position.
	
  cout << "+++ Test source position was defined as (RA,DEC) : ( ";
  if(fTestPosition)
    {
      cout << fTestRA << " , " << fTestDec << " ) " << endl;
      cout << "+++ RBM will be tuned to this test position. "
	   << endl;
      pfRingBackg = new VARingBackgroundAnl(fObsMode, fTestRA, fTestDec, pfOutputFile, fGroupID);
      pfHistogram2DAnl->setSourcePosition(fTestRA, fTestDec);
      pfWobbleAnl->
	setSourcePosition(fTestRA * TMath::RadToDeg(),
			  fTestDec * TMath::RadToDeg());
      if(!fExcludeCBGAnl)
	{
	  pfCrescentAnl->setSourcePosition(
					   fTestRA * TMath::RadToDeg(),
					   fTestDec * TMath::RadToDeg());
	  pfCrescentAnl->setSkyMap();
	}
		
    }
  else
    {
      cout << fSourceRA << " , " << fSourceDec << " ) " << endl;
      cout << "+++ RBM will be tuned to this test position. "
	   << endl;
      pfRingBackg = new VARingBackgroundAnl(fObsMode, fSourceRA, fSourceDec, pfOutputFile, fGroupID);
      pfHistogram2DAnl->setSourcePosition(fSourceRA, fSourceDec);
      pfWobbleAnl->
	setSourcePosition(fSourceRA * TMath::RadToDeg(),
			  fSourceDec * TMath::RadToDeg());
      if(!fExcludeCBGAnl)
	{
	  pfCrescentAnl->setSourcePosition(
					   fSourceRA * TMath::RadToDeg(),
					   fSourceDec * TMath::RadToDeg());
	}
    }
	
  //  fFileInput.close();
	
  //  fFileInput.open(FileListName.c_str());
	
  pair<double, double> position;
  double radius = fSourceExclusionRegionRadius;//[degree]
	
  if(fExcludeSource)
    {
	
      // set the exclusion region around the source.
		
      if(fTestPosition)
	{
	  position.first = fTestRA * 180. / TMath::Pi();
	  position.second = fTestDec * 180. / TMath::Pi();
	}
      else
	{
	  position.first = fSourceRA * 180. / TMath::Pi();
	  position.second = fSourceDec * 180. / TMath::Pi();
	}
      radius = fSourceExclusionRegionRadius; //degree
      // This is really really ugly but until I make sure position isn't used elsewhere it is the
      // easiest thing to do.
      VASkyMapExclusionRegion source((const char*)fSourceName.c_str(), position.first, position.second, radius);
		
      pfRingBackg->AddExclusionRegion(source);
      if(!fExcludeWobbleAnl)
	{
	  pfHistogram2DAnl->AddExclusionRegion(source);
	  pfWobbleAnl->AddExclusionRegion(source);
	}
      if(!fExcludeCBGAnl)
	{
	  pfCrescentAnl->AddExclusionRegion(source);
	}
    }
	
  if(fExcludeStars)
    {
      cout << "VAStage6Analysis: Generating list of star exclusion regions:" << endl;
      string catName = VAGlobal::instance()->getVEGASRuntimeDir() + "/common/share/starCat_mag8.root";
      if(fStarExclusionCatalog.length() != 0)
	{
	  catName = fStarExclusionCatalog;
	}
		
      VAStarCat cat(catName.c_str());
      if(!cat.isOK())
	{
	  cerr << "Yikes!  Couldn't open or load the catalog from " << catName << endl;
	  exit(-1);
	}
		
      //get things into degrees and grab a subcat around the pointing
      double s_ra   = fSourceRA * TMath::RadToDeg();
      double s_dec  = fSourceDec * TMath::RadToDeg();
      double radius = fStarCatFilterRadius;	  //if <= 0, no filtering is done
      double maglim = fStarMagLimit;
		
      // Star catalog should be centered at the test position if specified.  AJW--20090605
      if(fTestPosition)
	{
	  s_ra = fTestRA * TMath::RadToDeg();
	  s_dec = fTestDec * TMath::RadToDeg();
	}
		
      fStarCat = cat.getSubCatalog(s_ra, s_dec, radius, maglim);
      vector<VAStar> fovStars = fStarCat.getStars();
		
      cout << "VAStage6Analysis: Found " << fovStars.size() << " stars bright enough";
      cout << " to exclude within " << radius << " degrees of the source" << endl;
		
      for(size_t n = 0; n < fovStars.size(); n++)
	{
	  double xRad = fovStars.at(n).getExclusionRadius();
	  if(fStarExclusionRegionRadius > 0)
	    {
	      xRad = fStarExclusionRegionRadius;
	    }
			
	  VASkyMapExclusionRegion starExcl(fovStars.at(n).getBestName(),
					   fovStars.at(n).getRA(),
					   fovStars.at(n).getDec(),
					   xRad);
	  pfRingBackg->AddExclusionRegion(starExcl);
	  if(!fExcludeWobbleAnl)
	    {
	      pfHistogram2DAnl->AddExclusionRegion(starExcl);
	      pfWobbleAnl->AddExclusionRegion(starExcl);
	    }
	  if(!fExcludeCBGAnl)
	    {
	      pfCrescentAnl->AddExclusionRegion(starExcl);
	    }
	}
    }//if(fExcludeStars)
	
  if(fUserExclusionRegionsFile.length() > 0)
    {
	
      // Set the exclusion region for a userRegion.
      // Read the catalog of bright userRegions from the external file.
      // ------------------------------------------------------------------
      cout << "VAStage6Analysis: Reading user-defined exclusion regions from "
	   << fUserExclusionRegionsFile << endl;
      ifstream userRegionList;
      userRegionList.open(fUserExclusionRegionsFile.c_str());
      if(!userRegionList)
	{
	  // userRegion list file couldn't be opened
	  fErrstr.str("");
	  fErrstr << "Error: file " << fUserExclusionRegionsFile << " could not be opened!" << endl;
	  fErrstr << "+++ You won't be able to exclude these userRegions in analysis... " << endl;
	  VAGlobal::instance()->getErrorReport()
	    ->addError(fErrstr.str().c_str());
	}
      else
	{
		
	  double userRegionRA;
	  double userRegionDec;
	  double userRegionExclusionRadius;
	  char userRegionExclusionRegionNameChar[50], ch;
	  pair<double, double> userRegionPosition;
			
	  int i = 0;
	  while(!userRegionList.eof())
	    {
	      i++;
	      userRegionList >> userRegionRA
			     >> userRegionDec
			     >> userRegionExclusionRadius ;
	      //	 >> userRegionExclusionRegionName;
				
	      int j = 0;
	      userRegionExclusionRegionNameChar[0] = '\0';
	      while(userRegionList.get(ch))
		{
		  if(ch != '\n')
		    {
		      userRegionExclusionRegionNameChar[j++] = ch;
		    }
		  else
		    {
		      break;
		    }
		}
				
	      userRegionExclusionRegionNameChar[j++] = '\0';
				
	      string userRegionExclusionRegionName(userRegionExclusionRegionNameChar);
				
	      // 	  cout << "+++ " << i
	      // 	       << "--- Name : " << userRegionExclusionRegionName
	      // 	       << " r : " << userRegionExclusionRadius
	      // 	       << " RA : " << userRegionRA
	      // 	       << " Dec : " << userRegionDec << endl;
				
	      userRegionPosition.first = userRegionRA;
	      userRegionPosition.second = userRegionDec;
				
	      VASkyMapExclusionRegion userRegionCatalog(userRegionExclusionRegionName.c_str(),
							userRegionRA, userRegionDec, userRegionExclusionRadius);
						
	      pfRingBackg->AddExclusionRegion(userRegionCatalog);
	      if(!fExcludeWobbleAnl)
		{
		  pfHistogram2DAnl->AddExclusionRegion(userRegionCatalog);
		  pfWobbleAnl->AddExclusionRegion(userRegionCatalog);
		}
	      if(!fExcludeCBGAnl)
		{
		  pfCrescentAnl->AddExclusionRegion(userRegionCatalog);
		}
	    }
			
	  cout << "+++ Total of " << i << " exclusion regions has been read out! "
	       << endl;
	}
		
    } // end if ( fUseExclusionRegions )
	
  return true;
} //////////////////////////////////////end setup


bool VAStage6Analysis::performMainEventLoop()
{


  VARegisterThisFunction vreg(__PRETTY_FUNCTION__);
	
  int nFiles = 0;
  // major loop over the data files.
	
  TTree* pShowerTree = NULL;
  VARunHeader* pRunHeader = NULL;
  VAShowerData* pShowerData = new VAShowerData();
  VAParameterisedEventData* pParameterisedEventData = NULL;
	
  if(fReadFromStage5Combined)
    {
      // New objects for  size cutting (need access to Hillas data in combined tree)
      pParameterisedEventData = new VAParameterisedEventData();
    }
  double fBDTScore = 0;
	
  VAArrayInfo* pArrayInfo = NULL;
  VAQStatsData* pQStatsData = NULL;
  VAPixelStatusData* pPixelStatusData = NULL;
  double theta2 = 0; // to pass theta2 value from wobble analysis into
  // VAEventStatsTree.
  string option;
  float value;
  stringstream ss;
  long pPosition;
  string tmps;
  stringstream sstmps;
  uint32_t pGroupCounter = 0; // keeps track of the number of groups that have been processed
  uint32_t pGroupMemberCounter = 0; // number of run files processed in current group
  string RunFileName;
  string EAFileName;
	
  //fRunlists stores the contents of the runlist file
  fRunlists.clear();
  fRunlists.seekg(0, fRunlists.beg);
  std::size_t found;
  bool flag_ringSizeChange = false;
  bool runlistEndFlag = false;
  // start the infinite loop over the runlist file
  while(1)
    {
      // alternate between ON and OFF for On/Off obsMode
      bool runtype = true; //true = ON Run ; false = OFF Run
      if(fObsMode == "On/Off")
	{
	  if(nFiles % 2 == 0)
	    {
	      fRun = "On";
	    }
	  else
	    {
	      fRun = "Off";
	      runtype = false;
	    }
	}  // end of if(fObsMode=="On/Off")
		
      else if(fObsMode == "Wobble")
	{
	  fRun = "On";
	} //end if(fObsMode == "Wobble")
      ////////////////end On/Off obsMode
		
      getline(fRunlists, tmps); // read the first line
		
      if(tmps.substr(0, 9) == "[/RUNLIST")
	{
	  if(pGroupCounter == fNgroup) // if you are on the end of the last group
	    {
	      break; // from while(1) loop
	    }
	  getline(fRunlists, tmps);
	} //end if(tmps.substr(0,9) == "[/RUNLIST")
      if((tmps.substr(0, 1) == "[") && (tmps.substr(0, 2) != "[/"))
	{
	  //if there is only 1 run group, there is no need to process any further once a flag is seen. the EA and Config are still read
	  if(fNgroup == 1)
	    {
	      break; // from loop over runlistfile
	    }
	  //otherwise, continue and find the beginning of the next runlist group
			
	  sstmps.str(std::string());
	  sstmps.clear();
	  sstmps << "[RUNLIST ID: " << fGroupID.at(pGroupCounter) ;
	  found = fRunlists.str().find(sstmps.str());
	  cout << " finding runlist : " << sstmps.str() << "\t" << found << endl;
	  if(found < fRunlists.str().size() && fRunlists.good() && !fRunlists.eof())
	    {
	      fRunlists.seekg(fRunlists.beg + found);
	      getline(fRunlists, tmps); //reads [RUNLIST ID: line
	      getline(fRunlists, tmps); //reads the next runlist file
	      pGroupMemberCounter = 0; //initialization whenever reading new RUNLIST group
	      flag_ringSizeChange = false;
	    }
	  else
	    {
	      cout << "Cannot find runlist for group : " << fGroupID.at(pGroupCounter)
		   << " for " << pGroupCounter << endl;
	      break;
	    }
	} //end if((tmps.substr(0, 1) == "[") && (tmps.substr(0, 2) != "[/"))
		
      RunFileName = tmps;
		
      //nahee : very primitive way of finding right EA & configuration
      pPosition = fRunlists.tellg(); //bookmark position in file
      if(pGroupMemberCounter == 0) //if this is the first file in the group
	{
	  //find configuration
	  sstmps.str(std::string());
	  sstmps.clear();
	  sstmps << "[CONFIG ID: " << fGroupID.at(pGroupCounter) ;
	  found = fRunlists.str().find(sstmps.str());
	  if(found < fRunlists.str().size() && fRunlists.good() && !fRunlists.eof())
	    {
	      fRunlists.seekg(fRunlists.beg + found);
	      getline(fRunlists, tmps);
	      ss.str(std::string());
	      ss.clear();
	      pfStereoCutter->setToDefaultCuts();
	      while(true)
		{
		  // begin to loop over config parameters
		  getline(fRunlists, tmps);
		  ss << tmps;
		  if(tmps.substr(0, 8) == "[/CONFIG")
		    {
		      break; //from loop over config parameters
		    }
		  ss >> option >> value ; //read in the parameter name and value
		  ss.clear();
		  if(pfStereoCutter->setCutFromString(option, value) != 0)
		    {
		      if(option == "S6A_RingSize")
			{
			  if(value > 0 && value < 1.7)
			    {
			      fRingSize = value;
			      flag_ringSizeChange = true;
			    }
			}
		      else
			{
			  pfRingBackg->setCutsFromString(option, value, pGroupCounter);
			}
		    } //end if setting cuts from StereoCuts doesn't work
		  if(!fRunlists.good() && fRunlists.eof())
		    {
		      cout << "Footer search for CONFIG failed. Will break here" << endl;
		      break;
		    }
		} // end while loop over config lines
	    } // end if [CONFIG line is found
	  else if(pGroupCounter == 0)
	    {
	      //use command-line config for the first group.
	      cout << "Cannot find configuration for the first group. Will use command line configuration" << endl;
	    }
	  else
	    {
	      cout << "Cannot find configuration for group : " << fGroupID.at(pGroupCounter) << endl;
	      exit(-1);
	    }
			
			
	  if(fSpectrum || fUpperLimit)
	    {
	      if(fNgroup == 1 && fLookupFileName.size() > 0 && fLookupFileName.substr(fLookupFileName.size() - 5, 5) == ".root")
		{
		  EAFileName = fLookupFileName;
		}
	      sstmps.str(std::string()); //make sstmps blank
	      sstmps.clear(); // clear error flags, setting sstmps as goodbit
	      sstmps << "[EA ID: " << fGroupID.at(pGroupCounter);
	      found = fRunlists.str().find(sstmps.str()); //find EA ID flag
	      if(found < fRunlists.str().size() && fRunlists.good() && !fRunlists.eof())
		{
		  fRunlists.seekg(fRunlists.beg + found); //navigate to line signalling EA ID for current groupID
		  getline(fRunlists, tmps); //reads [EA ID: #] line
		  getline(fRunlists, tmps); //read the EA filename
		  if(tmps.substr(tmps.size() - 5, 5) == ".root")
		    {
		      EAFileName = tmps;
		      if(fNgroup == 1 && fLookupFileName.size() == 0 && fLookupFileName != EAFileName)
			{
			  cout << "++ You gave two EAs! We will use " << EAFileName
			       << " instead of " << fLookupFileName << endl;
			}
		    }
		  else
		    {
		      cout << "Error :: EA is not root file? Check that file extension is .root" << endl;
		      exit(-1);
		    }
		} // end if found EA configuration file
	      else
		{
		  //must find EA for any group except group 0
		  if((pGroupCounter == 0 && fLookupFileName == "") || pGroupCounter > 0)
		    {
		      cout << "Error :: Cannot find EA for Group ID : " << fGroupID.at(pGroupCounter) << endl;
		      cout << pGroupCounter << " " << fLookupFileName << " " << sDefaultLookupFileName << endl;
		      exit(-1);
		    }
		}
	    } //end if(fSpectrum || fUpperLimit)
	  pGroupCounter ++;
	} // end if(pGroupMemberCounter == 0)
      fRunlists.seekg(fRunlists.beg + pPosition); //pPosition bookmarks position after last processed file
      cout << " end chk :" << fRunlists.good() << " " << fRunlists.eof() << " " << fRunlists.tellg() << endl;
      if(!fRunlists.good() || fRunlists.eof() || runlistEndFlag)
	{
	  break; //break from main loop if the file is done
	}
      //       fFileInput >> RunFileName >> EAFileName;
      //       fFileInput.ignore(numeric_limits<streamsize>::max(),'\n');
      //       if (!fFileInput.good()) break;
      pGroupMemberCounter ++;
      std::cout << "\033[34m+++ Group ID : " << fGroupID.at(pGroupCounter - 1) << "  file # : " << pGroupMemberCounter << std::endl;
      std::cout << "\033[34m+++ " << nFiles + 1 << " Run file: " << RunFileName << "\033[m" << std::endl;
      if(fSpectrum || fUpperLimit)
	{
	  std::cout << "\033[34m+++ " << nFiles + 1 << " EA file: " << EAFileName << "\033[m" << std::endl;
	}
		
      if(pGroupMemberCounter == 1)
	{
	  std::cout << "\033[34m+++ Group Stereo Cut \033[m" << std::endl;
	  pfStereoCutter->printCuts();
	}
      if(flag_ringSizeChange && pGroupMemberCounter == 1)
	{
	  std::cout << "\033[34m+++ Setting Ring size to " << fRingSize <<" for group ID : " << fGroupID.at(pGroupCounter - 1) << std::endl;
	}
		
      // open the data file.
      VARootIO* pfInputRootIO = new VARootIO(RunFileName, true);
      //VAGlobal::instance()->addRootFile(&pfInputRootIO);
      pfInputRootIO->loadTheRootFile();
		
      if(!pfInputRootIO->IsOpen())
	{
	  throw VAFatalException("Couldn't open file " + RunFileName, __FILE__, __LINE__);
	}
		
      // open the EA file.
      VARootIO* pfEAInputRootIO = new VARootIO(EAFileName, true);
      VAEffectiveAreaManager* pfEAManager  = NULL;
      //VAGlobal::instance()->addRootFile(&pfEAInputRootIO);
      if(fSpectrum || fUpperLimit)
	{
	  pfEAManager = new VAEffectiveAreaManager();
	  pfEAInputRootIO->loadTheRootFile();
	  if(!pfEAInputRootIO->IsOpen())
	    {
	      throw VAFatalException("Couldn't open file " + EAFileName, __FILE__, __LINE__);
	    }
	  pfEAManager->loadEffectiveAreas(pfEAInputRootIO);
	  pfEAManager->getConfig().getValue("WindowForNoise", fWindowSizeForNoise);
	  cout << "Stage 6 Analysis picked up a window size for noise of " << fWindowSizeForNoise << endl;
			
			
	}
		
      if(fReadFromStage5Combined)
	{
	  pShowerTree = ((TTree*)pfInputRootIO->loadAnObject(gCombinedEventsTreeName, gSelectedEventsDirName, true));
	}
      else
	{
	  // Read events from stage4 or stage5
	  pShowerTree = (fReadFromStage4) ? ((TTree*)pfInputRootIO->loadAnObject(gShowerEventsTreeName, gShowerEventsDirName, true))
	    : ((TTree*)pfInputRootIO->loadAnObject(gShowerEventsTreeName, gSelectedEventsDirName, true));
	}
      if(!pShowerTree)
	throw VAFatalException("Couldn't load a shower tree from " + RunFileName + "\n"
			       "If you ran stage 5 on this file, then your file could be corrupt."
			       "If you did not run stage5 on this file then you can try reading the stage 4\n "
			       "data directly using the option -S6A_ReadFromStage4=true", __FILE__, __LINE__);
								   
      pShowerTree->SetBranchAddress(gShowerEventsBranchName.c_str(), &pShowerData);
      // Now need to add the Hillas guys in for the combined tree
      if(fReadFromStage5Combined)
	{
		
	  pShowerTree->SetBranchAddress(gCombinedEventsBDTScoreBranchName.c_str(), &fBDTScore);
	  pShowerTree->SetBranchAddress(gParEventsDataBranchName.c_str(), &pParameterisedEventData);
	}
		
      // Read info from run header.
      pRunHeader = pfInputRootIO->loadTheRunHeader();
      if(!pRunHeader)
	{
	  throw VAFatalException("Couldn't load the run header from " + RunFileName, __FILE__, __LINE__);
	}
		
      fRunNumber = pRunHeader->getRunNumber();
      fSourceID = pRunHeader->fRunInfo.fSourceID;
		
      double runStartMJD = pRunHeader->getStartTime().getMJDDbl();
      double runEndMJD = pRunHeader->getEndTime().getMJDDbl();
		
      std::cout << "+++ Run Number: " << fRunNumber << std::endl;
      std::cout << "+++ Source    : " << pRunHeader->fSourceInfo.fDescription.c_str() << std::endl;
      std::cout << "+++ RA        : " << pRunHeader->fSourceInfo.fRA << std::endl;
      std::cout << "+++ DEC       : " << pRunHeader->fSourceInfo.fDec << std::endl;
      std::cout << "+++ Epoch     : " << pRunHeader->fSourceInfo.fEpoch << std::endl;
      std::cout << "+++ RA Offset : " << pRunHeader->fRunInfo.fOffsetRA*(180. / TMath::Pi()) << std::endl;
      std::cout << "+++ DEC Offset: " << pRunHeader->fRunInfo.fOffsetDec*(180. / TMath::Pi()) << std::endl;
      std::cout << "+++ Source ID : " << pRunHeader->fRunInfo.fSourceID << std::endl;
      std::cout << "+++ Weather   : " << pRunHeader->fRunInfo.fWeather << std::endl;
		
      pArrayInfo = const_cast<VAArrayInfo*>(pfInputRootIO->loadTheArrayInfo(false));
      if(!pArrayInfo)
	{
	  throw VAFatalException("Couldn't load the array info from " + RunFileName, __FILE__, __LINE__);
	}
		
      VATime fStartTime = pRunHeader->getStartTime();
      VATime fEndTime = pRunHeader->getEndTime();
      int64_t exposureNS = fEndTime - fStartTime;
      double exposure = double(exposureNS) / 60. / 1.e9; // exposure = obs time from run header
      double fLiveTime = pRunHeader->pfRunDetails->fRunCutLiveTimeSeconds; // live time == effective live time!!!
		
      std::cout << "+++ Run Start Time : " << fStartTime << std::endl;
      std::cout << "+++ Run End Time   : " << fEndTime << std::endl;
      std::cout << "+++ Run exposure   : " << exposure  << " min" << std::endl;
      std::cout << "+++ Run Live Time  : " << fLiveTime / 60.  << " min" << std::endl;
		
      // insane sanity check on live time.
      if(fLiveTime < 1.)
	{
	  fLiveTime = exposure * 0.93 * 60.;
	  VAGlobal::instance()->getErrorReport()->addError("+++ Live Time is not defined, Run exposure will "
							   "be used instead and a 7% deadtime assumed.");
	}
		
      // increment the exposure for On or Off.
      if(fObsMode == "On/Off")
	{
	  if(nFiles % 2 == 0)
	    {
	      fExposureOn += fLiveTime / 60.;
	    }
	  else
	    {
	      fExposureOff += fLiveTime / 60.;
	    }
	}
      else
	{
	  fExposureOn += fLiveTime / 60.;
	}
		
      double sourceRA2000_Rad = 0;
      double sourceDec2000_Rad = 0;
      if(!fTestPosition)
	{
	  sourceRA2000_Rad = fSourceRA;
	  sourceDec2000_Rad = fSourceDec;
	}
      else
	{
	  sourceRA2000_Rad = fTestRA;
	  sourceDec2000_Rad = fTestDec;
	}
		
      // stereo analysis.
      double alpha = 1.0;
      double alphaCres = 1.0;
      double zenith = 0., minElevation = 90., maxZenith = 0., azimuthAtMaxZenith = 0;
      VATime timeAtMaxZenith;
      unsigned numEventsPassingCutsThisRun = 0;
      double wobbleOffsetRA2000_Rad  = 0.;
      double wobbleOffsetDec2000_Rad = 0.;
      double wobbleOffsetRA2000_Deg  = 0.;
      double wobbleOffsetDec2000_Deg = 0.;
      double wobbleOffset_Deg = 0.;
		
      if(pfShowerParamDisplay)
	{
	  int fNumEventEntries = (int)pShowerTree->GetEntries();
	  std::cout << "+++ ShowerData exists !     " << std::endl;
	  std::cout << "+++ STEREO: Number of events in the run : " << RunFileName << " : " << fNumEventEntries << std::endl;
	  pShowerTree->GetEntry(0);
			
	  // define the structure of stereo event.
	  pfDebug->DebugMessage(E_FULLDEBUG, "Getting shower data.");
	  // compute average zenith/azimuth and tracking directions.
	  pfDebug->DebugMessage(E_FULLDEBUG, "Computing average zenith.");
	  double trackingRA2000_Rad = 0, trackingDec2000_Rad = 0;
	  int counter = 0;
			
	  // ++++++ STARTING LOOP OVER EVENTS IN SINGLE RUN ++++++
	  for(int nEvent = 0; nEvent < fNumEventEntries; nEvent++)
	    {
	      pShowerTree->GetEntry(nEvent);
				
	      //write the event times for barycentering
	      if(fPulsarAnlFlag)
		{
		  VATime fArrayEventTime = pShowerData->fTime;
		  pfPulsar->fillBaryPhase(&fArrayEventTime, fSourceRA, fSourceDec);
		  //writes file for input for the TEMPO timing package
		  pfPulsar->writeTempoFile(&fArrayEventTime);
		}
				
	      if(pShowerData->fIsReconstructed)
		{
		  ++counter;
		  zenith += (90 - pShowerData->fArrayTrackingElevation_Deg);
		  if(pShowerData->fArrayTrackingElevation_Deg < minElevation)
		    {
		      minElevation = pShowerData->fArrayTrackingElevation_Deg;
		      azimuthAtMaxZenith = pShowerData->fArrayTrackingAzimuth_Deg;
		      timeAtMaxZenith = pShowerData->fTime;
		    }
		  trackingRA2000_Rad += pShowerData->fArrayTrackingRA_J2000_Rad;
		  trackingDec2000_Rad += pShowerData->fArrayTrackingDec_J2000_Rad;
		} // end if(pShowerData->fIsReconstructed)
	    } //++++++ END LOOP OVER EVENTS IN SINGLE RUN ++++++
			
			
	  //Barycenter and phasefold
	  if(fPulsarAnlFlag)
	    {
	      pfPulsar->barycenterWithTempo2();
	      pfPulsar->DoPhases();
	    }
			
	  maxZenith = 90 - minElevation;
	  zenith = zenith / (counter + 1.e-10);
	  trackingRA2000_Rad /= (counter + 1.e-10);
	  trackingDec2000_Rad /= (counter + 1.e-10);
			
	  std::cout << "+++ Ave Zenith angle : " << zenith << " [deg]" << std::endl;
	  std::cout << "+++ Max Zenith angle : " << maxZenith << " [deg]" << std::endl;
	  std::cout << "+++ Azimuth at Max Zenith : " << azimuthAtMaxZenith << " [deg]" << std::endl;
	  std::cout << "+++ Time at Max Zenith : " << timeAtMaxZenith.getString() << std::endl;
			
	  pfDebug->DebugMessage(E_FULLDEBUG, "Initializing the run for the Ring Background Model.");
	  if(!fExcludeRBMAnl)
	    {
	      pfRingBackg->ResetObservation(pRunHeader, trackingRA2000_Rad, trackingDec2000_Rad, pGroupCounter - 1);
	    }
	  pfDebug->DebugMessage(E_FULLDEBUG, "Initializing Run Stats");
	  if(!fExcludeRunStats) pfRunStats->reset(fRunNumber, fSourceID, runStartMJD, runEndMJD, (double)exposureNS / 1.e9,
						  fLiveTime, runtype, pArrayInfo);
	  pfDebug->DebugMessage(E_FULLDEBUG, "Determining test position and offset distance for the wobble analysis.");
	  wobbleOffsetRA2000_Rad = trackingRA2000_Rad - sourceRA2000_Rad;
	  wobbleOffsetDec2000_Rad = trackingDec2000_Rad - sourceDec2000_Rad;
	  wobbleOffset_Deg = slaSep(trackingRA2000_Rad, trackingDec2000_Rad, sourceRA2000_Rad, sourceDec2000_Rad) * TMath::RadToDeg();
	  wobbleOffsetRA2000_Deg = wobbleOffsetRA2000_Rad * TMath::RadToDeg();
	  wobbleOffsetDec2000_Deg = wobbleOffsetDec2000_Rad * TMath::RadToDeg();
			
	  std::cout << "+++ RA : " << sourceRA2000_Rad << " DEC : " << sourceDec2000_Rad << std::endl;
	  std::cout << "+++ WobbleAnl SourceRA     : " << sourceRA2000_Rad* TMath::RadToDeg() << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl SourceDec    : " << sourceDec2000_Rad* TMath::RadToDeg() << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl trackRA      : " << trackingRA2000_Rad* TMath::RadToDeg() << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl trackDec     : " << trackingDec2000_Rad* TMath::RadToDeg() << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl OffsetRA     : " << wobbleOffsetRA2000_Deg << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl OffsetDec    : " << wobbleOffsetDec2000_Deg << " deg" << std::endl;
	  std::cout << "+++ WobbleAnl Offset       : " << wobbleOffset_Deg << " deg" << std::endl;
			
	  pfDebug->DebugMessage(E_FULLDEBUG, "Initializing the wobble analysis.");
			
	  if(!fExcludeWobbleAnl)
	    {
	      pfWobbleAnl->setId(pRunHeader->getRunNumber(), pRunHeader->fRunInfo.fSourceID);
	      pfWobbleAnl->setNumbers(fNumRings, fRingSize);
	      pfWobbleAnl->setOffSet(wobbleOffsetRA2000_Deg, wobbleOffsetDec2000_Deg);
	      alpha = pfWobbleAnl->setBackgroundRegions();
	      pfHistogram2DAnl->setCenter(sourceRA2000_Rad, sourceDec2000_Rad);
	      pfHistogram2DAnl->setNumbers(fNumRings, fRingSize);
	      pfHistogram2DAnl->setBackgroundRegions(wobbleOffsetRA2000_Deg, wobbleOffsetDec2000_Deg);
	    }
	  if(!fExcludeCBGAnl)
	    {
	      pfCrescentAnl->setId(pRunHeader->getRunNumber(), pRunHeader->fRunInfo.fSourceID);
	      pfCrescentAnl->setOffsetPosition(wobbleOffsetRA2000_Deg, wobbleOffsetDec2000_Deg);
	      pfCrescentAnl->setTrackingPosition(trackingRA2000_Rad * TMath::RadToDeg(), trackingDec2000_Rad * TMath::RadToDeg());
	      pfCrescentAnl->setSkyMap();
	      pfCrescentAnl->setCrescent(fCBGRingLower, fCBGRingUpper);
	    }
			
	  // initialize run for effective areas
	  pQStatsData = const_cast<VAQStatsData*>(pfInputRootIO->loadTheQStatsData());
	  if(!pQStatsData)
	    {
	      throw VAFatalException("Couldn't load the qstats data from file " + RunFileName, __FILE__, __LINE__);
	    }
			
	  pPixelStatusData = pfInputRootIO->loadThePixelStatusData();
	  if(!pPixelStatusData)
	    {
	      throw VAFatalException("Couldn't load the pixel status data!", __FILE__, __LINE__);
	    }
			
	  // initialize run for spectrum / UL.
	  if(fSpectrum || fUpperLimit)
	    {
	      pfDebug->DebugMessage(E_FULLDEBUG, "Initializing Spectrum");
				
	      // initialize the EA parameter data
	      pfEAParameterData->pfQStatsData = pQStatsData;
	      pfEAParameterData->pfPixelStatusData = pPixelStatusData;
	      pfEAParameterData->pfRunHeader = pRunHeader;
	      pfEAParameterData->pfArrayInfo = pArrayInfo;
	      pfEAParameterData->pfShowerData = pShowerData;
	      pfEAParameterData->fAzimuth = azimuthAtMaxZenith;
	      pfEAParameterData->fZenith = maxZenith;
	      pfEAParameterData->fTime = timeAtMaxZenith;
	      fLTParams = pfEAManager->getParametersForThisEvent(pfEAParameterData);
				
	      double Threshold = 0;
	      double EnergyMax = 0;
	      // Still need to implement this function...
	      // Remeber that this is before we loop over the events... So we'll need to use a time-independent estimator
	      // We need tp use:
	      //  - maxZenith
	      //  - azimuthAtMaxZenith
	      //  - offset
	      //  - noise...
	      pfDebug->DebugMessage(E_FULLDEBUG, "Azimuth angle (deg) used to get safe energy range: ", azimuthAtMaxZenith);
	      pfDebug->DebugMessage(E_FULLDEBUG, "Zenith angle (deg) used to get safe energy range: ", maxZenith);
	      // Now call the super function.
	      float minSafeEnergy_GeV = 0;
	      float maxSafeEnergy_GeV = 0;
	      if(pfEAManager->getSafeEnergyRange(fLTParams, wobbleOffset_Deg, minSafeEnergy_GeV, maxSafeEnergy_GeV))
		{
		  pfDebug->DebugMessage(E_FULLDEBUG, "Determined safe energy threshold: ", minSafeEnergy_GeV);
		}
	      else
		{
		  pfDebug->DebugMessage(E_FULLDEBUG, "Failed to determine safe energy range, you should be worried!");
		}
				
	      Threshold = (double)minSafeEnergy_GeV / 1000.0;
	      EnergyMax = (double)maxSafeEnergy_GeV / 1000.0;
	      pfRunStats->fillSpectrumStats((double)minSafeEnergy_GeV,
					    (double)maxSafeEnergy_GeV);
	      pfSpectrum->InitNewRun(fObsMode, fRun, Threshold, EnergyMax);
	      pfLiveTime->InitNewRun(fObsMode, fRun, Threshold, EnergyMax);
	    } // end of if(fSpectrum / UL)
	  // ++++ LOOP OVER EVENTS +++++ AGAIN!?!?!?!?!?!?!
	  pfDebug->DebugMessage(E_FULLDEBUG, "Starting loop over events.");
	  for(int nEvent = 0; nEvent < fNumEventEntries; nEvent++)
	    {
	      pShowerTree->GetEntry(nEvent);
	      // Modify the energy scale if required
	      if(fApplyEnergyCorrection)
		{
		  TFormula tempFormula("tempFormula", fEnergyCorrectionFormula.c_str());
		  float energy_temp = (float)tempFormula.Eval(pShowerData->fEnergy_GeV);
		  // same RMS as percentage value.
		  // is it correct? maybe not, but maybe yes. NG
		  pShowerData->fEnergyRMS_GeV *= energy_temp / pShowerData->fEnergy_GeV;
		  pShowerData->fEnergy_GeV = energy_temp;
		}
	      vector<int> fTelIndex;
	      fTelIndex.clear();
	      // Fill RunStats for raw events (no cuts).
	      if(!fExcludeRunStats)
		{
		  pfRunStats->fillRaw(pShowerData, fRun);
		}
	      // apply simple cuts for mean scaled Width/Length.
	      bool passSizeCut = true;
				
	      if(fReadFromStage5Combined && pShowerData->fIsReconstructed)
		{
		  passSizeCut = !pfSizeMaskCutter->isEventCutBySizeMask(pParameterisedEventData);
		}
	      if((fBDTScore >= -1 && passSizeCut && pShowerData->fIsReconstructed
		  && !(pfStereoCutter->isCut(pShowerData, fCutMask, fBDTScore) && fCutMask != VA_FAILED_THETA)
		  && !pfTelMaskCutter->isEventCutByReconTelescopeMask(pShowerData)) ||
		 (fBDTScore < -1 && passSizeCut && pShowerData->fIsReconstructed
		  && !(pfStereoCutter->isCut(pShowerData, fCutMask) && fCutMask != VA_FAILED_THETA)
		  && !pfTelMaskCutter->isEventCutByReconTelescopeMask(pShowerData)))
		{
		  pfDebug->DebugMessage(E_FULLDEBUG, "Event passes all cuts");
		  ++numEventsPassingCutsThisRun;
		  // Initialization of the flag.
		  string flag = "NULL";
		  string flagw = "NULL";
		  string flagC = "NULL";
		  double eventWeightC = 1.;
		  double eventWeight = 1.;
					
		  // in the "Ring Background Model" the background region is defined
		  // as the region similar to on region but reflected with respect to
		  // the center of the camera.
		  if(!fExcludeRBMAnl)
		    {
		      pfRingBackg->fill(pShowerData, fRun, flag);
		      pfDebug->DebugMessage(E_FULLDEBUG, "Flag from RBM:  ", flag);
		    }
		  // Fill RunStats class for cut events.
		  if(!fExcludeRunStats)
		    {
		      pfRunStats->fillCut(pShowerData, fRun, flag);
		    }
		  if(fMCData)
		    {
		      flag = "On";
		    }
					
		  // here is the call to the 'display' class. the flag comes from the RBM.
		  if(fDisplayStereoParamFlag && flag != "NULL")
		    {
		      pfShowerParamDisplay->fillShower(fMCData, pArrayInfo, pShowerData, flag, fTelIndex);
		    }
					
		  // event selection in the Wobble mode.
		  if(fObsMode == "Wobble")
		    {
		      string fRunW = "Wobble";
						
						
		      // This is the method for the On/Off event selection.
		      // flagw is just the same flag but for wobble.
		      if(!fExcludeWobbleAnl)
			{
			  pfHistogram2DAnl->getEventOrBg(pShowerData, wobbleOffsetRA2000_Deg, wobbleOffsetDec2000_Deg, flagw, eventWeight, theta2);
			  pfHistogram2DAnl->fill(pShowerData, fRunW);
			  pfWobbleAnl->fillEvents(pShowerData);
			}
						
		      if(!fExcludeCBGAnl)
			{
			  pfCrescentAnl->setEventPosition(pShowerData->fDirectionRA_J2000_Rad * TMath::RadToDeg(), pShowerData->fDirectionDec_J2000_Rad * TMath::RadToDeg());
			  flagC = pfCrescentAnl->getEventOrBg();
			  eventWeightC = pfCrescentAnl->fillEvents(flagC);
			  pfCrescentAnl->fillHist();
			}
						
		    } // end of if(fObsMode == "Wobble")
		  // now define the flag and weight.
		  if(fBgEstimate == "wobble")
		    {
		      flag = flagw; // Switch flag, eventWeight comes from Wobble.
		      if(flag == "on")
			{
			  flag = "On";
			}
		      if(flag == "off")
			{
			  flag = "Off";
			}
		    }
		  else if(fBgEstimate == "crescent")
		    {
		      flag = flagC;
		      eventWeight = eventWeightC;
		    }
		  else if(fBgEstimate == "ring")
		    {
		      eventWeight = 1.;   // Use RBM flag, redefine eventWeight.
		    }
		  pfDebug->DebugMessage(E_FULLDEBUG, "Flag for this event:  ", flag);
		  if(flag == "On")
		    {
		      eventWeight = 1.;    // Always single signal region.
		    }
		  pfDebug->DebugMessage(E_FULLDEBUG, "Weight for this event:  ", eventWeight);
					
					
					
		  // if in signal region (thetasqr) select the event for the pulsar analysis
		  // and redefine on/off on=event has phase that falls into pulsed region
		  // off=event has phase that falls into phase region used to estimate the background
		  //  else do not use the event in the spectrum reconstructuion =NULL
		  if(fPulsarAnlFlag) //need redefinve eventWeight
		    {
		      if(flag == "On")
			{
			  // select event for pulsar analysis
			  pfPulsar->SelectEventForAnalysis(nEvent,  pShowerData->fEnergy_GeV);
			  if(pfPulsar->IsOn(nEvent))
			    {
			      flag = "On";
			    }
			  else if(pfPulsar->IsOff(nEvent))
			    {
			      flag = "Off";
			    }
			  else
			    {
			      flag = "NULL";
			    }
			}
		      else
			{
			  flag = "NULL";
			}
		    } //end if(fPulsarAnlFlag)
					
					
		  // TO BE FIXED!!! Temporary solution for scope of instructions. NG
		  //pfEAParameterData->fAzimuth= pShowerData->fArrayTrackingAzimuth_Deg;
		  //pfEAParameterData->fZenith=90.-pShowerData->fArrayTrackingElevation_Deg;
		  pfEAParameterData->fAzimuth = (float)TMath::RadToDeg() * pShowerData->fDirectionAzimuth_Rad;
		  pfEAParameterData->fZenith = 90.0 - (float)TMath::RadToDeg() * pShowerData->fDirectionElevation_Rad;
		  pfEAParameterData->fTime = pShowerData->fTime;
		  if(fSpectrum || fUpperLimit)
		    {
		      fLTParams = pfEAManager->getParametersForThisEvent(pfEAParameterData);
		    }
					
		  float Area = 0;
		  // fill the SPECTRUM histograms.
		  if(fSpectrum && flag != "NULL")
		    {
		      //double sigmaarea = 0;
		      pfEAManager->fill_MigrationMatrix(pfEAParameterData);
						
						
		      if(fNumOutOfBoundsForArea.first.size() != fLTParams.size())
			{
			  for(uint32_t dim = 0; dim < fLTParams.size(); dim++)
			    {
			      fNumOutOfBoundsForArea.first.push_back(0);
			    }
			}
		      if(fNumOutOfBoundsForArea.second.size() != fLTParams.size())
			{
			  for(uint32_t dim = 0; dim < fLTParams.size(); dim++)
			    {
			      fNumOutOfBoundsForArea.second.push_back(0);
			    }
			}
						
		      float Error = 0;
		      if(pfEAManager->getInterpolatedEffectiveAreaAndErrorFromTable(fLTParams,
										    TMath::Log10(pShowerData->fEnergy_GeV / 1000), 0, Area, Error, pfEAParameterData,
										    fNumOutOfBoundsForArea.first, fNumOutOfBoundsForArea.second))
			{
			  pfSpectrum->fillEvent(pShowerData, wobbleOffset_Deg, (double)Area, (double)Error, flag);
			  pfSpectrum->FillEffectiveArea(pfEAManager->getEffectiveAreaCurve(fLTParams));
			}
		    } // end of if( fSpectrum ) block
					
		  // fill the UPPER LIMIT object, this is quite a mess-around due to the bad design of stage6, an exclusive OR on Wobble and RBM must be implemented.
		  if(fUpperLimit && flag != "NULL")
		    {
		      if(fBgEstimate == "ring" && flag == "On")
			{
			  pfUpperLimit->FillEffectiveArea(pfEAManager->getEffectiveAreaCurve(fLTParams));
			}
		      else if(fBgEstimate == "wobble" || fBgEstimate == "crescent")
			{
			  pfUpperLimit->FillEvent(pShowerData, flag);
			  pfUpperLimit->FillEffectiveArea(pfEAManager->getEffectiveAreaCurve(fLTParams));
			}
		    }
					
					
		  if(!fExcludeEventStats)
		    {
					
		      double phase = 0.0;
		      double MJDBary = 0.0;
		      if(fPulsarAnlFlag)
			{
			  phase = pfPulsar->GetPhase(nEvent);
			  MJDBary =  pfPulsar->GetMJDBary(nEvent);
			}
						
		      bool bOff = (flag == "Off") ? true : false;
		      bool bOn  = (flag == "On") ? true : false;
						
		      pfEventStats->fill(pShowerData, bOn, bOff, eventWeight, getAverageNoise(pRunHeader, pQStatsData, pShowerData->fTime), Area, theta2, phase, MJDBary, pfEAParameterData->fOffset, fBDTScore);
						
		    } //end if(!fExcludeEventStats)
		} // pfStereoCutter and pTelMaskCutter
	    } // end of the event loop.
	  pfDebug->DebugMessage(E_FULLDEBUG, "Finished event loop");
	  pfDebug->DebugMessage(E_FULLDEBUG, "Closed event dumps");
	} // end of if ShowerData exists and stereo analysis enabled
		
		
      // Add the migration matrix
      if(fSpectrum)
	{
	  if(pfMigrationMatrix->GetEntries() == 0)
	    {
	      Int_t nBinsX = pfEAManager->getMigrationMatrix()->GetNbinsX();
	      const double* xBins = pfEAManager->getMigrationMatrix()->GetXaxis()->GetXbins()->GetArray();
	      Int_t nBinsY = pfEAManager->getMigrationMatrix()->GetNbinsY();
	      const double* yBins = pfEAManager->getMigrationMatrix()->GetYaxis()->GetXbins()->GetArray();
	      pfMigrationMatrix->SetBins(nBinsX, xBins, nBinsY, yBins);
	      std::cout << "+++ Migration Matrix histogram binning set +++" << endl;
	    }
	  pfMigrationMatrix->Add(pfEAManager->getMigrationMatrix());
	}
		
      // Estimate relative exposure for this run.
      pfDebug->DebugMessage(E_FULLDEBUG, "Getting runwise acc curve");
      double relExposure = 0;
      TH1F* accplot = 0;
      if(!fExcludeRBMAnl)
	{
	  accplot = pfRingBackg->endOfRun();
	}
      pfDebug->DebugMessage(E_FULLDEBUG, "Estimating relative exposure");
      if(fDoRelativeExposure)
	{
	  relExposure = pfExposureCalc->calculateExposure(zenith, wobbleOffset_Deg, fRingSize, accplot);
	}
      else
	//relExposure = numEventsPassingCutsThisRun;
	{
	  relExposure = 1;
	}
		
      std::cout << "Run " << fRunNumber << " Relative Exposure = " << relExposure << std::endl;
		
      // wobble analysis after a single run.
      pfDebug->DebugMessage(E_FULLDEBUG, "post-run wobble analysis");
      if(fObsMode == "Wobble")
	{
	  if(!fExcludeWobbleAnl)
	    {
	      pfWobbleAnl->setHist(getFilledHistogram(pfHistogram2DAnl));
	      //	pfWobbleAnl->draw();
	      std::cout << "+++ Wobble RA offset : " << wobbleOffsetRA2000_Deg << " DEC offset : "  << wobbleOffsetDec2000_Deg << std::endl;
	      pfDebug->DebugMessage(E_FULLDEBUG, "wobble increment exposure");
	      pfWobbleAnl->addEventsPerRun(fLiveTime / 60., relExposure);
	      // BJZ hack:
	      alpha = pfWobbleAnl->getAlpha(0);
	      //cout << "alpha: " << alpha << " zn Correction factor: " << pfWobbleAnl->getZnCorrectionFactor() << endl;
	      pfWobbleAnl->write(pfOutputFile);
	      pfDebug->DebugMessage(E_FULLDEBUG, "fill2d increment exposure");
	      pfHistogram2DAnl->incrementExposure(relExposure);
	    }
	  if(!fExcludeCBGAnl)
	    {
	      alphaCres = pfCrescentAnl->getAlpha();
	      pfCrescentAnl->addEventsPerRun(fLiveTime / 60., relExposure);
	    }
			
	}
		
      // spectrum analysis after a single run.
      pfDebug->DebugMessage(E_FULLDEBUG, "post-run spectral analysis");
      float alphaPfSpectrum = 1.0; // appropriate for 'RBM' which has 1 BG region
      if(fBgEstimate == "wobble")
	{
	  alphaPfSpectrum = alpha;    // if using Wobble for BG instead
	}
      if(fBgEstimate == "crescent")
	{
	  alphaPfSpectrum = alphaCres;    // if using CBG for BG instead
	}
      if(fPulsarAnlFlag)
	{
	  alphaPfSpectrum = pfPulsar->GetSignalBackgroundRatio();
	}
      pfLiveTime->incrementExposure(fObsMode, fRun, fLiveTime / 60.*60., relExposure, alphaPfSpectrum);
      pfDebug->DebugMessage(E_FULLDEBUG, "Finished all stereo analysis for this run");
      // end of stereo analysis.
		
      // Write RunStats tree entries for this run.
      if(!fExcludeRunStats)
	{
	  pfRunStats->fillTree();
	}
      nFiles++;
		
      pPixelStatusData->Delete();
      pPixelStatusData = NULL;
      pQStatsData->Delete();
      pQStatsData = NULL;
      pArrayInfo->Delete();
      pArrayInfo = NULL;
      pRunHeader->Delete();
      pRunHeader = NULL;
      pShowerTree->Delete();
      pShowerTree = NULL;
      if(fSpectrum)
	{
	  cout << endl
	       << "----------------------------------------------" << endl
	       << "Report from effective area out-of-bounds check" << endl
	       << "----------------------------------------------" << endl
	       << "Number of events that fell outside the range covered by the"
	       << " effective area file, by dimension AND parameter. If a"
	       << " numberis greater than zero, the effective area file may not"
	       << " be appropriate for all data in your analysis. Consider"
	       << " using a file with better coverage in that dimension."
	       << endl;
	  cout << endl;
	  cout << "Dimension: area (numbers given as below/above range)" << endl;
			
	  bool caughtOutOfBounds = false;
			
	  for(uint32_t dim = 0; dim < pfEAManager->fParams.size(); dim++)
	    {
	      std::string dimName = pfEAManager->VAEffectiveAreaDimensionTypeToString(pfEAManager->fParams.at(dim).first);
	      cout << dimName;
	      cout << " " << fNumOutOfBoundsForArea.first.at(dim) << "/" << fNumOutOfBoundsForArea.second.at(dim);
	      cout << endl;
				
	      if(pfEAManager->fParams.at(dim).first != E_AbsoluteOffset)
		{
		  if((fNumOutOfBoundsForArea.first.at(dim) > 0) ||
		     (fNumOutOfBoundsForArea.second.at(dim) > 0))
		    {
		      caughtOutOfBounds = true;
		    }
		}
	    } // end for loop over dim
	  cout << endl;
			
	  if(caughtOutOfBounds)
	    {
	      ostringstream errstr;
	      errstr.str("");
	      errstr << "Warning! At least one event fell outside the range"
		     << " covered by the effective area file. Use the out-of-"
		     << "bounds report output to guide your selection of the"
		     << " appropriate file.";
	      VAGlobal::instance()->getErrorReport()->addError(errstr.str().c_str());
	    }
			
	} //end if(fSpectrum)
		
      if(pfEAManager)
	{
	  delete pfEAManager; // Need to delete before pfInputRootIO is closed
	  pfEAManager = NULL;
	}
		
      if(pfInputRootIO)
	{
	  pfInputRootIO->closeTheRootFile();
	  delete pfInputRootIO;
	  pfInputRootIO = NULL;
	}
		
      if(pfEAInputRootIO)
	{
	  pfEAInputRootIO->closeTheRootFile();
	  delete pfEAInputRootIO;
	  pfEAInputRootIO = NULL;
	}
    } // end of the loop over runlist file
  // --
  // stereo analysis results.
	
  std::cout << std::endl << std::endl << std::endl << "Finished loop over files - final steps:" << std::endl << std::endl << std::endl;
	
  pfDebug->DebugMessage(E_FULLDEBUG, "Post loop over files pulsar analysis");
  if(fPulsarAnlFlag)
    {
      cout << endl << "Here starts the final pulsar analysis, good luck" << endl;
      cout << "================================================" << endl << endl;
		
      if(fSpectrum)
	{
	  pfPulsar->SetEnergyBins(pfSpectrum->GetRebinnedSpectrumHist());
	}
      pfPulsar->CreatePhaseogramAndCalculateSignificance();
      cout << endl;
      pfPulsar->writePhaseogram(pfOutputFile);
      pfPulsar->writeFileWithTimesOfSelectedEvents();
      cout << endl << "That was everything from the pulsar analysis, hope you are happy with the results" << endl << endl;
    } //end if(fPulsarAnlFlag)
	
	
  if(pfShowerParamDisplay)
    {
	
      // Write EventStats tree to file
      pfDebug->DebugMessage(E_FULLDEBUG, "write eventstats");
      if(!fExcludeEventStats)
	{
	  pfEventStats->writeToFile(pfOutputFile);
	}
		
      pfDebug->DebugMessage(E_FULLDEBUG, "write wobble");
      if(fObsMode == "Wobble")
	{
	  if(!fExcludeWobbleAnl)
	    {
	      pfHistogram2DAnl->draw(fBatchMode);
	      pfHistogram2DAnl->write();
	    }
	  if(!fExcludeCBGAnl)
	    {
	      pfCrescentAnl->draw();
	    }
	}
      // final analysis for the RBM.
      pfDebug->DebugMessage(E_FULLDEBUG, "rbm last steps");
      // 	if(fObsMode=="Wobble" && !fExcludeRBMAnl ){
      if(!fExcludeRBMAnl)
	{
	  pfRingBackg->ConfigureAcceptanceMap(fSaveAcceptance, fLoadAcceptance,
					      fConfigDirectory, fAcceptPlotLookup);
	  pfRingBackg->MakeAcceptanceMap(fBatchMode);
	}
		
      pfDebug->DebugMessage(E_FULLDEBUG, "rbm write");
      if(!fExcludeRBMAnl && !fExcludeRBMFinal)
	{
	  pfRingBackg->MakeSignificanceMap(fExposureOn, fBatchMode);
	}
      pfRingBackg->write();
		
      if(fDisplayStereoParamFlag && (!fBatchMode))
	{
	  pfShowerParamDisplay->draw();
	}
      if(fDisplayStereoParamFlag)
	{
	  pfShowerParamDisplay->write(pfOutputFile);
	}
      TDirectory* d = (TDirectory*)pfOutputFile->mkdir("LiveTime");
      d->cd();
      pfLiveTime->Write("VALiveTime");
		
      // final analysis for the spectrum evaluation.
      pfDebug->DebugMessage(E_FULLDEBUG, "write spectrum");
      char* wrdirname = (char*)"IndecentResults";
      if(fSpectrum)
	{
	  pfSpectrum->Finish();
	  if(!fBatchMode)
	    {
	      pfSpectrum->Draw();
	    }
	  char* dirname = (fApplyEnergyCorrection) ? (wrdirname) : ((char*)"Spectrum");
	  if(!pfOutputFile)
	    {
	      cout << "ERROR: unable to write output file" << endl;
	    }
	  if(!pfOutputFile->Get(dirname))
	    {
	      d = (TDirectory*)pfOutputFile->mkdir(dirname);
	    }
	  d->cd();
	  char* spname = (fApplyEnergyCorrection) ? ((char*)"IndecentSpectrum") : ((char*)"VASpectrumAnl");
	  pfSpectrum->Write(spname);
	  char* ufname = (fApplyEnergyCorrection) ? ((char*)"IndecentUnfolding") : ((char*)"Unfolding");
	  if(!!pfOutputFile->Get(ufname))
	    {
	      d = (TDirectory*)pfOutputFile->mkdir(ufname);
	    }
	  d->cd();
	  pfMigrationMatrix->Scale(1. / pfMigrationMatrix->Integral());
	  pfMigrationMatrix->Write("hMigrationMatrix");
	  pfSpectrum->GetFullExcessHist()->Write("hFullExcessHist");
	}
		
      if(fUpperLimit)
	{
	  char* dirname = (fApplyEnergyCorrection) ? (wrdirname) : ((char*)"UpperLimit");
	  if(!pfOutputFile->Get(dirname))
	    {
	      d = (TDirectory*)pfOutputFile->mkdir(dirname);
	    }
	  d->cd();
	  if(fBgEstimate == "ring")
	    {
	      pfUpperLimit->CalculateUpperLimit((int)pfRingBackg->getNon(), (int)pfRingBackg->getNoff(), pfRingBackg->getRBMAlpha());
	    }
	  else if(fBgEstimate == "wobble" || fBgEstimate == "crescent")
	    {
	      pfUpperLimit->CalculateUpperLimit();
	    }
	  char* ulname = (fApplyEnergyCorrection) ? ((char*)"IndecentUpperLimit") : ((char*)"VAUpperLimit");
	  pfUpperLimit->Write(ulname);
	  if(!fBatchMode)
	    {
	      pfUpperLimit->Draw();
	    }
	}
      // draw RBM plots.
      pfDebug->DebugMessage(E_FULLDEBUG, "draw rbm");
      if(!fBatchMode)
	{
	  if(!fExcludeRBMAnl && !fExcludeRBMFinal)
	    {
	      pfRingBackg->draw(fObsMode, fStarCat, fDrawExclusionRegions);
				
	      ////////////////////////////////////////////////////////////
	      // if(fDrawExclusionRegions > 1 && fStarCat.isOK()){
	      //               vector<TCanvas*> *plots = pfRingBackg->getPlots();
	      //               for(unsigned n=0;n<plots->size();n++){
	      //                 if(plots) fStarCat.draw((TPad*)plots->at(n), , true);
	      //               }
	      //             } else {
	      //               cerr << "Starcat not loaded!  Cannot overplot stars!" << endl;
	      //             }
	      pfRingBackg->print(fObsMode, fOutputFileName, fConfigDirectory);
	    }
	}//if(!fBatchMode){
    }//if(pfShowerParamDisplay){
	
  if(fStarCat.isOK())
    {
      //Can write the thing to the output file, too if you want
      if((TDirectory*)pfOutputFile->cd("RingBackgroundModelAnalysis"))
	{
	  fStarCat.Write("fStarCat");
	}
    }
	
  // Write the Runwise Statistics tree.
  pfDebug->DebugMessage(E_FULLDEBUG, "write runstats");
  if(!fExcludeRunStats)
    {
      pfRunStats->write(pfOutputFile);
    }
  pfDebug->DebugMessage(E_FULLDEBUG, "done file writing");
	
	
  if(fLightCurve)
    {
      /*
	TGraphErrors *lightCurve = pfLightCurveAnl->plotLightCurve(
	if(!fBatchMode)lightCurve->Draw("AP");
      */
    }
	
  //
  // Output out-of-bounds report to the log
  //
	
  if(fObsMode == "Wobble")
    {
      if(!fExcludeWobbleAnl)
	{
	  pfWobbleAnl->calcSigni();
	}
      if(!fExcludeCBGAnl)
	{
	  pfCrescentAnl->calcSignificance();
	}
    }
	
  VAGlobal::instance()->getErrorReport()->printErrorReport();
	
  cout << "+++ Stage 6 : resultExtractor has finished! " << endl;
  if(fBatchMode)
    {
      exit(EXIT_SUCCESS);
    }
  cout << "+++ In order to quit the program please click 'File' button "
       << endl;
  cout << "+++ in any of pop-up windows on the screen and then 'Quit ROOT'"
       << endl;
  return true;
} // end performMainEventLoop function

Double_t VAStage6Analysis::getAverageNoise(VARunHeader* pRunHeader, VAQStatsData* pQStatsData, const VATime& time)
{
  if(!pRunHeader)
    {
      throw VAFatalException("The run header is null", __FILE__, __LINE__);
    }
	
  if(!pRunHeader->pfRunDetails)
    {
      throw VAFatalException("The run details is null", __FILE__, __LINE__);
    }
	
  if(!pQStatsData)
    {
      throw VAFatalException("The qstats is null", __FILE__, __LINE__);
    }
	
  Double_t meanNoise = 0;
  Int_t nCameras = 0;
	
  if(fWindowSizeForNoise == 0)
    {
      cout << "Stage 6 Analysis picked up a window size of 0, which likely means no EA is available.  Reverting to default size 7 for noise calc in EventStatsTree.\n";
      fWindowSizeForNoise = 7;
    }
  for(int telID = 0; telID < kMaxTels; telID++)
    {
      if((pRunHeader->pfRunDetails->fExpectedTels.at(telID)) == 1)
	{
	  float noise = pQStatsData->getCameraAverageTraceVar(telID, fWindowSizeForNoise, time);
	  if(noise > 0)
	    {
	      meanNoise += noise;
	      nCameras++;
	    }
	}
    }// end for loop over telID
	
  if(nCameras > 0)
    {
      meanNoise /= nCameras;
    }
	
  return meanNoise;
} //end getAverageNoise function

