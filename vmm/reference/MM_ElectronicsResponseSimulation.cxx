/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 *  MM_ElectronicsResponseSimulation.cxx
 *  MC for micromegas athena integration
 *
 **/
#include "MM_Digitization/MM_ElectronicsResponseSimulation.h"

#include <numeric>

#include "AthenaKernel/getMessageSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "TF1.h"

#include "CLHEP/Random/RandGaussZiggurat.h"

/*******************************************************************************/
MM_ElectronicsResponseSimulation::MM_ElectronicsResponseSimulation(ConfigModule&& module):
    m_cfg{std::move(module)} {
    // let the VMM peak search run in a wider window in order to simulate the dead time and to avoid a bug at the upper limit of the time
    // window
    m_vmmShaper =
        std::make_unique<VMM_Shaper>(m_cfg.peakTime, m_cfg.timeWindowLowerOffset - m_cfg.vmmDeadtime, m_cfg.timeWindowUpperOffset + m_cfg.vmmUpperGrazeWindow);
    m_vmmShaper->initialize();
}
/*******************************************************************************/
MM_DigitToolOutput MM_ElectronicsResponseSimulation::getPeakResponseFrom(const MM_ElectronicsToolInput& digiInput) const {
    DataCache cache{};
    vmmPeakResponseFunction(cache, digiInput);

    /// ToDo: include loop for calculating Trigger study vars
    // MM_DigitToolOutput(bool hitWasEff, std::vector <int> strpos, std::vector<float> time, std::vector<int> charge, std::vector<float>
    // threshold, int strTrig, float strTimeTrig ):
    MM_DigitToolOutput to_ret{true, std::move(cache.nStripElectronics), 
                                 std::move(cache.tStripElectronicsAbThr), 
                                 std::move(cache.qStripElectronics), 5, 0.3};

    return to_ret;
}
/*******************************************************************************/
MM_DigitToolOutput MM_ElectronicsResponseSimulation::getThresholdResponseFrom(const MM_ElectronicsToolInput& digiInput) const {
    DataCache cache{};
    vmmThresholdResponseFunction(cache, digiInput);
    MM_DigitToolOutput to_ret{true, std::move(cache.nStripElectronics), 
                                std::move(cache.tStripElectronicsAbThr), 
                                std::move(cache.qStripElectronics), 5, 0.3};
    return to_ret;
}
/*******************************************************************************/
void MM_ElectronicsResponseSimulation::vmmPeakResponseFunction(DataCache& cache, const MM_ElectronicsToolInput& digiInput) const {
    
    
    const std::vector<int>& numberofStrip = digiInput.NumberOfStripsPos();
    const std::vector<std::vector<float>>& qStrip =  digiInput.chipCharge();
    const std::vector<std::vector<float>>& tStrip = digiInput.chipTime();
    const std::vector<float>& stripsElectronicsThreshold = digiInput.stripThreshold();
    
    for (unsigned int ii = 0; ii < numberofStrip.size(); ii++) {
        // find min and max times for each strip:
        bool thisStripFired{false};
        double leftStripFired = false;
        double rightStripFired = false;

        // find the maximum charge:
        if (m_cfg.useNeighborLogic) {  // only check neighbor strips if VMM neighbor logic is enabled
            if (ii > 0) {
                leftStripFired =
                    m_vmmShaper->hasChargeAboveThreshold(qStrip.at(ii - 1), tStrip.at(ii - 1), stripsElectronicsThreshold.at(ii - 1));
            }

            if (ii + 1 < numberofStrip.size()) {
                rightStripFired =
                    m_vmmShaper->hasChargeAboveThreshold(qStrip.at(ii + 1), tStrip.at(ii + 1), stripsElectronicsThreshold.at(ii + 1));
            }
        }

        thisStripFired = m_vmmShaper->hasChargeAboveThreshold(qStrip.at(ii), tStrip.at(ii), stripsElectronicsThreshold.at(ii));

        // check if neighbor strip was above threshold
        bool neighborFired = leftStripFired || rightStripFired;

        // Look at strip if it or its neighbor was above threshold  and if neighbor logic of the VMM is enabled:
        if (thisStripFired || (m_cfg.useNeighborLogic && neighborFired)) {
            double charge{FLT_MAX}, time{FLT_MAX};

            // if strip is below threshold but read through NL reduce threshold to low value of 1e to still find the peak
            float tmpScaledThreshold = (thisStripFired ? stripsElectronicsThreshold.at(ii) : 1);

            bool foundPeak = m_vmmShaper->vmmPeakResponse(qStrip.at(ii), tStrip.at(ii), tmpScaledThreshold, charge, time);
            if (!foundPeak) continue;  // if no peak was found within the enlarged time window the strip is skipped and no digit is created.
            if (time < m_cfg.timeWindowLowerOffset || time > m_cfg.timeWindowUpperOffset)
                continue;  // only accept strips in the correct time window

            cache.nStripElectronics.push_back(numberofStrip.at(ii));
            cache.tStripElectronicsAbThr.push_back(time);
            cache.qStripElectronics.push_back(charge);
        }
    }
}

void MM_ElectronicsResponseSimulation::vmmThresholdResponseFunction(DataCache& cache, const MM_ElectronicsToolInput& digiInput) const {
    
    const std::vector<int>& numberofStrip = digiInput.NumberOfStripsPos();
    const std::vector<std::vector<float>>& qStrip = digiInput.chipCharge();
    const std::vector<std::vector<float>>& tStrip = digiInput.chipTime();
    const std::vector<float>& electronicsThreshold = digiInput.stripThreshold();

    for (unsigned int ii = 0; ii < numberofStrip.size(); ii++) {
        double localThresholdt{FLT_MAX}, localThresholdq{FLT_MAX};

        bool crossedThreshold =
            m_vmmShaper->vmmThresholdResponse(qStrip.at(ii), tStrip.at(ii), electronicsThreshold.at(ii), localThresholdq, localThresholdt);
        if (!crossedThreshold)
            continue;  // if no threshold crossing was found within the time window the strip is skipped and no digits are created.
        if (localThresholdt < m_cfg.timeWindowLowerOffset || localThresholdt > m_cfg.timeWindowUpperOffset)
            continue;  // only accept strips in the correct time window

        cache.nStripElectronics.push_back(numberofStrip.at(ii));
        cache.tStripElectronicsAbThr.push_back(localThresholdt);
        cache.qStripElectronics.push_back(localThresholdq);
    }
}
