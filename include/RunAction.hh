#pragma once
// RunAction.hh
// Handles ROOT TTree output (with CSV fallback) per worker thread.
// Each worker writes its own file; collect_results.py runs hadd + analysis.

#include "G4UserRunAction.hh"
#include "SimConfig.hh"
#include "EventData.hh"

#include <string>
#include <memory>

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction(const SimConfig& cfg, bool isMaster);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    // Called by EventAction on each worker thread
    void RecordEvent(const EventData& data);

private:
    // One-row-per-worker provenance tree; see include/RunMeta.hh for why it
    // has to live in the file rather than in the submit script's manifest.
    // Written at end of run, because `thrown` is not known before then.
    void WriteRunMeta();

    const SimConfig& fConfig;
    bool             fIsMaster;

    struct Impl;
    std::unique_ptr<Impl> fImpl;

    // Run accumulators (for summary printout)
    long   fTotalEvents   = 0;
    double fSumEdepDrift  = 0.0;
    double fSumNprimDrift = 0.0;
    double fSumEdepAmp    = 0.0;
    double fSumNprimAmp   = 0.0;
    long   fWrittenEvents = 0;   // rows actually written (--skip-empty)
    long   fSkippedEmpty  = 0;
};
