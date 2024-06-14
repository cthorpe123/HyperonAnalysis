#ifndef _EventAssembler2_h_
#define _EventAssembler2_h_

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "SimParticle.h"
#include "RecoParticle.h"
#include "Event.h"

using std::vector;
using std::string;

///////////////////////////////////////////////////////////////////////////////////////////////

class EventAssembler {

   private:

      // Input file and event tree
      TFile * f_in=nullptr;
      TTree * t_in;

      TTree *t_meta;
      Double_t POT;

      // Number of events (entries in event tree)
      int nEvents;

      TString DataDir;

   public:

      EventAssembler();
      EventAssembler(bool loadweights);

      // Setters and getters
      void SetFile(string infilename,string sampletype);
      void Close();

      Event GetEvent(int i);
      Long64_t GetNEvents();
      double GetPOT();		

   private:

      std::string SampleType;
      const bool LoadWeights;

      // General Info
      Bool_t          IsData;
      UInt_t          EventID;
      Int_t           run;
      Int_t           subrun;
      Int_t           event;
      Double_t        Weight;

      // Event flags
      vector<string>* Mode=0;
      vector<string>* CCNC=0;
      Int_t NMCTruths;
      Int_t NMCTruthsInTPC;
      vector<bool>* InActiveTPC=0;
      vector<bool>* IsHyperon=0;
      vector<bool>* IsLambda=0;
      vector<bool>* IsLambdaCharged=0;
      vector<bool>* IsSigmaZero=0;
      vector<bool>* IsSigmaZeroCharged=0;
      vector<bool>* IsAssociatedHyperon=0;
      vector<bool>* IsSignal=0;
      vector<bool>* IsSignalSigmaZero=0;
      Bool_t GoodReco;
      Bool_t EventHasNeutronScatter;
      Bool_t EventHasHyperon;

      // Truth variables
      std::vector<SimParticle>* Neutrino=0;
      std::vector<SimParticle>* Lepton=0;
      std::vector<SimParticle>* Hyperon=0;
      std::vector<SimParticle>* PrimaryNucleon=0;
      std::vector<SimParticle>* PrimaryPion=0;
      std::vector<SimParticle>* PrimaryKaon=0;
      std::vector<SimParticle>* Decay=0;
      std::vector<SimParticle>* SigmaZeroDecayLambda=0;
      std::vector<SimParticle>* SigmaZeroDecayPhoton=0;
      std::vector<SimParticle>* KaonDecay=0;

      vector<double>  *TruePrimaryVertex_X=0;
      vector<double>  *TruePrimaryVertex_Y=0;
      vector<double>  *TruePrimaryVertex_Z=0;
      vector<double>  *DecayVertex_X=0;
      vector<double>  *DecayVertex_Y=0;
      vector<double>  *DecayVertex_Z=0;

      // Reco variables
      TVector3* RecoPrimaryVertex=0;
      Int_t NPrimaryTrackDaughters;
      Int_t NPrimaryShowerDaughters;
      std::vector<RecoParticle>* TracklikePrimaryDaughters=0;
      std::vector<RecoParticle>* ShowerlikePrimaryDaughters=0;

      vector<vector<int>>* ConnSeedIndexes_Plane0=0;
      vector<vector<int>>* ConnOutputIndexes_Plane0=0;
      vector<vector<int>>* ConnOutputSizes_Plane0=0;
      vector<vector<int>>* ConnSeedChannels_Plane0=0;
      vector<vector<int>>* ConnSeedTicks_Plane0=0;
      vector<vector<int>>* ConnSeedIndexes_Plane1=0;
      vector<vector<int>>* ConnOutputIndexes_Plane1=0;
      vector<vector<int>>* ConnOutputSizes_Plane1=0;
      vector<vector<int>>* ConnSeedChannels_Plane1=0;
      vector<vector<int>>* ConnSeedTicks_Plane1=0;
      vector<vector<int>>* ConnSeedIndexes_Plane2=0;
      vector<vector<int>>* ConnOutputIndexes_Plane2=0;
      vector<vector<int>>* ConnOutputSizes_Plane2=0;
      vector<vector<int>>* ConnSeedChannels_Plane2=0;
      vector<vector<int>>* ConnSeedTicks_Plane2=0;

      // Systematics
      vector<string>* SysDials=0;
      vector<vector<vector<double>>>* SysWeights=0;

      std::vector<RecoParticle>* RepassTracklikePrimaryDaughters=0;
};

///////////////////////////////////////////////////////////////////////////////////////////////

#endif
