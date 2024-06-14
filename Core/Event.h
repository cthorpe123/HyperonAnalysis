#ifndef _Event2_h_
#define _Event2_h_

#include <vector>

#include "TVector3.h"

#include "SimParticle.h"
#include "RecoParticle.h"

using std::vector;
using std::string;

struct Event {

	bool IsData=false;
	UInt_t EventID;
	Int_t run;
	Int_t subrun;
	Int_t event;
	Int_t fileID;
	double Weight=1.0;

        // Flags applying to the entire event
        Bool_t EventIsSignal = false;
        Bool_t EventHasHyperon = false;
        Bool_t EventHasNeutronScatter = false;
        Bool_t EventIsSignalSigmaZero = false;

        // Flags for each MCTruth
        vector<string> Mode;
        vector<string> CCNC;
        Int_t NMCTruths;
        Int_t NMCTruthsInTPC;
        vector<bool> InActiveTPC;
        vector<bool> IsHyperon;
        vector<bool> IsLambda;
        vector<bool> IsLambdaCharged;
        vector<bool> IsSigmaZero;
        vector<bool> IsSigmaZeroCharged;
        vector<bool> IsAssociatedHyperon;
        vector<bool> IsSignal;
        vector<bool> IsSignalSigmaZero;
        Bool_t GoodReco;
                
        vector<TVector3> TruePrimaryVertex;

	//true variables
	std::vector<SimParticle> Neutrino;
	std::vector<SimParticle> Lepton;
	std::vector<SimParticle> Hyperon;
	std::vector<SimParticle> PrimaryNucleon;
	std::vector<SimParticle> PrimaryPion;
	std::vector<SimParticle> PrimaryKaon;
	std::vector<SimParticle> Decay;
	std::vector<SimParticle> SigmaZeroDecayLambda;
	std::vector<SimParticle> SigmaZeroDecayPhoton;
	std::vector<SimParticle> KaonDecay;

        vector<TVector3> DecayVertex;

	TVector3 RecoPrimaryVertex;
	Int_t NPrimaryTrackDaughters;
	Int_t NPrimaryShowerDaughters;

	std::vector<RecoParticle> TracklikePrimaryDaughters;
	std::vector<RecoParticle> ShowerlikePrimaryDaughters;

	int TrueMuonIndex = -1;
	int TrueDecayProtonIndex = -1;
	int TrueDecayPionIndex = -1;

	// Muon, proton and pion candidates
	RecoParticle MuonCandidate;
	RecoParticle DecayProtonCandidate;
	RecoParticle DecayPionCandidate;

	Float_t SelectorBDTScore;
	Float_t AnalysisBDTScore;

        vector<vector<int>> ConnSeedIndexes_Plane0;
        vector<vector<int>> ConnOutputIndexes_Plane0;
        vector<vector<int>> ConnOutputSizes_Plane0;
        vector<vector<int>> ConnSeedChannels_Plane0;
        vector<vector<int>> ConnSeedTicks_Plane0;
        vector<vector<int>> ConnSeedIndexes_Plane1;
        vector<vector<int>> ConnOutputIndexes_Plane1;
        vector<vector<int>> ConnOutputSizes_Plane1;
        vector<vector<int>> ConnSeedChannels_Plane1;
        vector<vector<int>> ConnSeedTicks_Plane1;
        vector<vector<int>> ConnSeedIndexes_Plane2;
        vector<vector<int>> ConnOutputIndexes_Plane2;
        vector<vector<int>> ConnOutputSizes_Plane2;
        vector<vector<int>> ConnSeedChannels_Plane2;
        vector<vector<int>> ConnSeedTicks_Plane2;
    
        vector<string> SysDials;
        vector<vector<vector<double>>> SysWeights;

	std::vector<RecoParticle> RepassTracklikePrimaryDaughters;

        void Print(){ std::cout << run << "  " << subrun << "  " << event << std::endl; }

};

#endif
