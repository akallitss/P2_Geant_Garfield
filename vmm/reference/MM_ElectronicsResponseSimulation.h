/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MM_DIGITIZATION_MM_ELECTRONICSRESPONSESIMULATION_H
#define MM_DIGITIZATION_MM_ELECTRONICSRESPONSESIMULATION_H

#include <TF1.h>

/// STD'S
#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

/// Projects
#include "MM_Digitization/MM_DigitToolOutput.h"
#include "MM_Digitization/MM_ElectronicsToolInput.h"

// VMM Mapping
#include "MM_Digitization/VMM_Shaper.h"
namespace CLHEP{
  class HepRandomEngine;
}
class MM_ElectronicsResponseSimulation {
public:
    
    struct ConfigModule {
      /** power of responce function */
      float peakTime{0.};
      float timeWindowLowerOffset{0.};
      float timeWindowUpperOffset{0.};
      float stripDeadTime{0.};
      float artDeadTime{0.};
      float vmmDeadtime{0.};
      float vmmUpperGrazeWindow{0.};
      bool useNeighborLogic{false};
    };

    MM_ElectronicsResponseSimulation(ConfigModule&& module);
    ~MM_ElectronicsResponseSimulation() = default;   
   

    MM_DigitToolOutput getPeakResponseFrom(const MM_ElectronicsToolInput& digiInput) const;
    MM_DigitToolOutput getThresholdResponseFrom(const MM_ElectronicsToolInput& digiInput) const;
   
    
    float getPeakTime() const { return m_cfg.peakTime; };
    float getTimeWindowLowerOffset() const { return m_cfg.timeWindowLowerOffset; };
    float getTimeWindowUpperOffset() const { return m_cfg.timeWindowUpperOffset; };
    float getVmmDeadTime() { return m_cfg.vmmDeadtime; }
    float getVmmUpperGrazeWindow() { return m_cfg.vmmUpperGrazeWindow; }
    float getStripdeadtime() const { return m_cfg.stripDeadTime; };
    float getARTdeadtime() const { return m_cfg.artDeadTime; };

  

private:
    const ConfigModule m_cfg{};
    
    struct DataCache {
        std::vector<float> tStripElectronicsAbThr{};
        std::vector<float> qStripElectronics{};
        std::vector<int> nStripElectronics{};
        std::unique_ptr<VMM_Shaper> vmmShaper{};
    };
    
    std::unique_ptr<VMM_Shaper> m_vmmShaper{};
      
    void vmmPeakResponseFunction(DataCache& cache, const MM_ElectronicsToolInput& digiInput) const;
    void vmmThresholdResponseFunction(DataCache& cache, const MM_ElectronicsToolInput& digiInput) const;

};

#endif
