#ifndef _EventAssembler_cxx_
#define _EventAssembler_cxx_

#include "EventAssembler.h"

EventAssembler::EventAssembler() :
LoadWeights(true)
{

}

EventAssembler::EventAssembler(bool loadweights) :
LoadWeights(loadweights)
{

}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load a file containing ntuples

void EventAssembler::SetFile(string infilename,string sampletype){

   SampleType = sampletype;

   DataDir = getenv("DATA_DIR");         

   TString file = DataDir + TString("/") +  TString(infilename);

   f_in = TFile::Open( file );

   f_in->GetObject("ana/OutputTree",t_in);
   //f_in->GetObject("HyperonNTuples/OutputTree",t_in);

   nEvents = t_in->GetEntries();

   t_in->SetBranchStatus("Neutrino",1);
   t_in->SetBranchStatus("Lepton",1);
   t_in->SetBranchStatus("Hyperon",1);
   t_in->SetBranchStatus("Decay",1);
   t_in->SetBranchStatus("PrimaryNucleon",1);
   t_in->SetBranchStatus("PrimaryPion",1);
   t_in->SetBranchStatus("PrimaryKaon",1);
   t_in->SetBranchStatus("SigmaZeroDecayLambda",1);
   t_in->SetBranchStatus("SigmaZeroDecayPhoton",1);
   t_in->SetBranchStatus("KaonDecay",1);
   t_in->SetBranchStatus("TracklikePrimaryDaughters",1);
   t_in->SetBranchStatus("ShowerlikePrimaryDaughters",1);
   t_in->SetBranchStatus("IsData",1);
   t_in->SetBranchStatus("EventID",1);
   t_in->SetBranchStatus("run",1);
   t_in->SetBranchStatus("subrun",1);
   t_in->SetBranchStatus("event",1);
   t_in->SetBranchStatus("Weight",1);
   t_in->SetBranchStatus("Mode",1);
   t_in->SetBranchStatus("CCNC",1);
   t_in->SetBranchStatus("NMCTruths",1);
   t_in->SetBranchStatus("NMCTruthsInTPC",1);
   t_in->SetBranchStatus("InActiveTPC",1);
   t_in->SetBranchStatus("IsHyperon",1);
   t_in->SetBranchStatus("IsLambda",1);
   t_in->SetBranchStatus("IsLambdaCharged",1);
   t_in->SetBranchStatus("IsSigmaZero",1);
   t_in->SetBranchStatus("IsSigmaZeroCharged",1);
   t_in->SetBranchStatus("IsAssociatedHyperon",1);
   t_in->SetBranchStatus("IsSignal",1);
   t_in->SetBranchStatus("IsSignalSigmaZero",1);
   t_in->SetBranchStatus("GoodReco",1);
   t_in->SetBranchStatus("EventHasNeutronScatter",1);
   t_in->SetBranchStatus("EventHasHyperon",1);
   t_in->SetBranchStatus("TruePrimaryVertex_X",1);
   t_in->SetBranchStatus("TruePrimaryVertex_Y",1);
   t_in->SetBranchStatus("TruePrimaryVertex_Z",1);
   t_in->SetBranchStatus("DecayVertex_X",1);
   t_in->SetBranchStatus("DecayVertex_Y",1);
   t_in->SetBranchStatus("DecayVertex_Z",1);
   t_in->SetBranchStatus("RecoPrimaryVertex",1);
   t_in->SetBranchStatus("NPrimaryTrackDaughters",1);
   t_in->SetBranchStatus("NPrimaryShowerDaughters",1);
   t_in->SetBranchStatus("ConnSeedIndexes_Plane0",1);
   t_in->SetBranchStatus("ConnOutputIndexes_Plane0",1);
   t_in->SetBranchStatus("ConnOutputSizes_Plane0",1);
   t_in->SetBranchStatus("ConnSeedChannels_Plane0",1);
   t_in->SetBranchStatus("ConnSeedTicks_Plane0",1);
   t_in->SetBranchStatus("ConnSeedIndexes_Plane1",1);
   t_in->SetBranchStatus("ConnOutputIndexes_Plane1",1);
   t_in->SetBranchStatus("ConnOutputSizes_Plane1",1);
   t_in->SetBranchStatus("ConnSeedChannels_Plane1",1);
   t_in->SetBranchStatus("ConnSeedTicks_Plane1",1);
   t_in->SetBranchStatus("ConnSeedIndexes_Plane2",1);
   t_in->SetBranchStatus("ConnOutputIndexes_Plane2",1);
   t_in->SetBranchStatus("ConnOutputSizes_Plane2",1);
   t_in->SetBranchStatus("ConnSeedChannels_Plane2",1);
   t_in->SetBranchStatus("ConnSeedTicks_Plane2",1);
   if(LoadWeights){
   t_in->SetBranchStatus("SysDials",1);
   t_in->SetBranchStatus("SysWeights",1);
   }
   t_in->SetBranchAddress("IsData", &IsData);
   t_in->SetBranchAddress("EventID", &EventID);
   t_in->SetBranchAddress("run", &run);
   t_in->SetBranchAddress("subrun", &subrun);
   t_in->SetBranchAddress("event", &event);
   t_in->SetBranchAddress("Weight", &Weight);
   t_in->SetBranchAddress("Mode", &Mode);
   t_in->SetBranchAddress("CCNC", &CCNC);
   t_in->SetBranchAddress("NMCTruths", &NMCTruths);
   t_in->SetBranchAddress("NMCTruthsInTPC", &NMCTruthsInTPC);
   t_in->SetBranchAddress("InActiveTPC", &InActiveTPC);
   t_in->SetBranchAddress("IsHyperon", &IsHyperon);
   t_in->SetBranchAddress("IsLambda", &IsLambda);
   t_in->SetBranchAddress("IsLambdaCharged", &IsLambdaCharged);
   t_in->SetBranchAddress("IsSigmaZero", &IsSigmaZero);
   t_in->SetBranchAddress("IsSigmaZeroCharged", &IsSigmaZeroCharged);
   t_in->SetBranchAddress("IsAssociatedHyperon", &IsAssociatedHyperon);
   t_in->SetBranchAddress("IsSignal", &IsSignal);
   t_in->SetBranchAddress("IsSignalSigmaZero", &IsSignalSigmaZero);
   t_in->SetBranchAddress("GoodReco", &GoodReco);
   t_in->SetBranchAddress("EventHasNeutronScatter", &EventHasNeutronScatter);
   t_in->SetBranchAddress("EventHasHyperon", &EventHasHyperon);

   t_in->SetBranchAddress("TruePrimaryVertex_X", &TruePrimaryVertex_X);
   t_in->SetBranchAddress("TruePrimaryVertex_Y", &TruePrimaryVertex_Y);
   t_in->SetBranchAddress("TruePrimaryVertex_Z", &TruePrimaryVertex_Z);

   t_in->SetBranchAddress("Neutrino",&Neutrino);
   t_in->SetBranchAddress("Lepton",&Lepton);
   t_in->SetBranchAddress("Hyperon",&Hyperon);
   t_in->SetBranchAddress("Decay",&Decay);
   t_in->SetBranchAddress("PrimaryNucleon",&PrimaryNucleon);
   t_in->SetBranchAddress("PrimaryPion",&PrimaryPion);
   t_in->SetBranchAddress("PrimaryKaon",&PrimaryKaon);
   t_in->SetBranchAddress("SigmaZeroDecayLambda",&SigmaZeroDecayLambda);
   t_in->SetBranchAddress("SigmaZeroDecayPhoton",&SigmaZeroDecayPhoton);
   t_in->SetBranchAddress("KaonDecay",&KaonDecay);

   t_in->SetBranchAddress("DecayVertex_X", &DecayVertex_X);
   t_in->SetBranchAddress("DecayVertex_Y", &DecayVertex_Y);
   t_in->SetBranchAddress("DecayVertex_Z", &DecayVertex_Z);

   t_in->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
   t_in->SetBranchAddress("NPrimaryTrackDaughters", &NPrimaryTrackDaughters);
   t_in->SetBranchAddress("NPrimaryShowerDaughters", &NPrimaryShowerDaughters);
   t_in->SetBranchAddress("TracklikePrimaryDaughters",&TracklikePrimaryDaughters);
   t_in->SetBranchAddress("ShowerlikePrimaryDaughters",&ShowerlikePrimaryDaughters);
               
   t_in->SetBranchAddress("ConnSeedIndexes_Plane0", &ConnSeedIndexes_Plane0);
   t_in->SetBranchAddress("ConnOutputIndexes_Plane0", &ConnOutputIndexes_Plane0);
   t_in->SetBranchAddress("ConnOutputSizes_Plane0", &ConnOutputSizes_Plane0);
   t_in->SetBranchAddress("ConnSeedChannels_Plane0", &ConnSeedChannels_Plane0);
   t_in->SetBranchAddress("ConnSeedTicks_Plane0", &ConnSeedTicks_Plane0);
   t_in->SetBranchAddress("ConnSeedIndexes_Plane1", &ConnSeedIndexes_Plane1);
   t_in->SetBranchAddress("ConnOutputIndexes_Plane1", &ConnOutputIndexes_Plane1);
   t_in->SetBranchAddress("ConnOutputSizes_Plane1", &ConnOutputSizes_Plane1);
   t_in->SetBranchAddress("ConnSeedChannels_Plane1", &ConnSeedChannels_Plane1);
   t_in->SetBranchAddress("ConnSeedTicks_Plane1", &ConnSeedTicks_Plane1);
   t_in->SetBranchAddress("ConnSeedIndexes_Plane2", &ConnSeedIndexes_Plane2);
   t_in->SetBranchAddress("ConnOutputIndexes_Plane2", &ConnOutputIndexes_Plane2);
   t_in->SetBranchAddress("ConnOutputSizes_Plane2", &ConnOutputSizes_Plane2);
   t_in->SetBranchAddress("ConnSeedChannels_Plane2", &ConnSeedChannels_Plane2);
   t_in->SetBranchAddress("ConnSeedTicks_Plane2", &ConnSeedTicks_Plane2);
   if(LoadWeights){
   t_in->SetBranchAddress("SysDials", &SysDials);
   t_in->SetBranchAddress("SysWeights", &SysWeights);
   }

   t_in->SetBranchAddress("RepassTracklikePrimaryDaughters",&RepassTracklikePrimaryDaughters);

   // Get the metadata tree
   f_in->GetObject("ana/MetaTree",t_meta);
   t_meta->SetBranchAddress("POT",&POT);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Close the file

void EventAssembler::Close(){

   if(f_in != nullptr) f_in->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the POT used to generate the current file

double EventAssembler::GetPOT(){

   int n_meta = t_meta->GetEntries();

   double TotalPOT=0.0;

   for(int i_meta=0;i_meta<n_meta;i_meta++){
      t_meta->GetEntry(i_meta);
      TotalPOT += POT;
   }

   return TotalPOT;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the number of events in the current file

Long64_t EventAssembler::GetNEvents(){

   return nEvents;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load an event

Event EventAssembler::GetEvent(int i){

   Event e;

   t_in->GetEntry(i);

   e.IsData = IsData;
   e.run = run;
   e.subrun = subrun;
   e.event = event;

   e.Weight = Weight;
   e.Mode = *Mode;

   // Special cases of "interaction mode"
   if(SampleType == "Data") e.Mode = {"Data"};
   else if(SampleType == "EXT") e.Mode = {"EXT"};
   else if(SampleType == "Dirt") e.Mode = {"Dirt"};

   e.CCNC = *CCNC;
   e.InActiveTPC = *InActiveTPC;
   e.IsHyperon = *IsHyperon;
   e.IsLambda = *IsLambda;
   e.IsLambdaCharged = *IsLambdaCharged;
   e.IsSigmaZero = *IsSigmaZero;
   e.IsSigmaZeroCharged = *IsSigmaZeroCharged;
   e.IsSignal = *IsSignal;
   e.IsSignalSigmaZero = *IsSignalSigmaZero;
   e.IsAssociatedHyperon = *IsAssociatedHyperon;
   e.GoodReco = GoodReco;
   e.EventHasNeutronScatter = EventHasNeutronScatter;
   e.EventHasHyperon = EventHasHyperon;

   e.NMCTruths = NMCTruths;
   e.NMCTruthsInTPC = NMCTruthsInTPC;

   e.TruePrimaryVertex.clear();

   for(size_t i_v=0;i_v<TruePrimaryVertex_X->size();i_v++)
      e.TruePrimaryVertex.push_back(TVector3(TruePrimaryVertex_X->at(i_v),TruePrimaryVertex_Y->at(i_v),TruePrimaryVertex_Z->at(i_v)));

   e.Neutrino = *Neutrino;
   e.Lepton = *Lepton;
   e.Hyperon = *Hyperon;
   e.PrimaryNucleon = *PrimaryNucleon;
   e.PrimaryPion = *PrimaryPion;
   e.PrimaryKaon = *PrimaryKaon;
   e.Decay = *Decay;
   e.SigmaZeroDecayLambda = *SigmaZeroDecayLambda;
   e.SigmaZeroDecayPhoton = *SigmaZeroDecayPhoton;
   e.KaonDecay = *KaonDecay;

   e.DecayVertex.clear();

   for(size_t i_v=0;i_v<DecayVertex_X->size();i_v++)
      e.DecayVertex.push_back(TVector3(DecayVertex_X->at(i_v),DecayVertex_Y->at(i_v),DecayVertex_Z->at(i_v)));

   e.RecoPrimaryVertex = *RecoPrimaryVertex;
   e.NPrimaryTrackDaughters = NPrimaryTrackDaughters;
   e.NPrimaryShowerDaughters = NPrimaryShowerDaughters;

   e.TracklikePrimaryDaughters = *TracklikePrimaryDaughters;
   e.ShowerlikePrimaryDaughters = *ShowerlikePrimaryDaughters;

   e.ConnSeedIndexes_Plane0 = *ConnSeedIndexes_Plane0;
   e.ConnOutputIndexes_Plane0 = *ConnOutputIndexes_Plane0;
   e.ConnOutputSizes_Plane0 = *ConnOutputSizes_Plane0;
   e.ConnSeedChannels_Plane0 = *ConnSeedChannels_Plane0;
   e.ConnSeedTicks_Plane0 = *ConnSeedTicks_Plane0;
   e.ConnSeedIndexes_Plane1 = *ConnSeedIndexes_Plane1;
   e.ConnOutputIndexes_Plane1 = *ConnOutputIndexes_Plane1;
   e.ConnOutputSizes_Plane1 = *ConnOutputSizes_Plane1;
   e.ConnSeedChannels_Plane1 = *ConnSeedChannels_Plane1;
   e.ConnSeedTicks_Plane1 = *ConnSeedTicks_Plane1;
   e.ConnSeedIndexes_Plane2 = *ConnSeedIndexes_Plane2;
   e.ConnOutputIndexes_Plane2 = *ConnOutputIndexes_Plane2;
   e.ConnOutputSizes_Plane2 = *ConnOutputSizes_Plane2;
   e.ConnSeedChannels_Plane2 = *ConnSeedChannels_Plane2;
   e.ConnSeedTicks_Plane2 = *ConnSeedTicks_Plane2;

   if(LoadWeights){
   e.SysDials = *SysDials;
   e.SysWeights = *SysWeights;
   }

   e.RepassTracklikePrimaryDaughters = *RepassTracklikePrimaryDaughters;

   return e;
}

///////////////////////////////////////////////////////////////////////////////////////////////

#endif
