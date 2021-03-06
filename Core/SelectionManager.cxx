#ifndef _SelectionManager_cxx_
#define _SelectionManager_cxx_

#include "SelectionManager.h"

///////////////////////////////////////////////////////////////////////////////////////////////

// Default constructor

SelectionManager::SelectionManager() : 
   a_FluxWeightCalc(0) ,
   a_GenG4WeightCalc() ,
   a_FiducialVolume(1,0.0) ,
   a_MuonID(-0.35,18.0,0.0) ,
   a_TrackLengthCutManager(1000.0,1000.0) , 
   a_SelectorBDTManager("Test") ,
   a_AnalysisBDTManager("Test") ,
   a_EventListFilter() ,
   a_SecondaryVertexFitter(5) , 
   a_CTTest_Plane0(0) ,
   a_CTTest_Plane1(1) ,
   a_CTTest_Plane2(2) 
{
   DeclareCuts();
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Destructor

SelectionManager::~SelectionManager(){


}

///////////////////////////////////////////////////////////////////////////////////////////////

// Constructor that loads selection parameters

SelectionManager::SelectionManager(SelectionParameters p) :
   a_FluxWeightCalc(p.p_BeamMode) ,
   a_GenG4WeightCalc() ,
   a_FiducialVolume(p.p_FV,p.p_Padding) ,
   a_MuonID(p.p_PID_Cut,p.p_Minimum_MIP_Length,p.p_Max_Displacement) ,
   a_TrackLengthCutManager(p.p_SecondaryTrackLengthCut,p.p_SecondaryTrackLengthCut) ,
   a_SelectorBDTManager("Test") , 
   a_AnalysisBDTManager("Test") ,
   a_EventListFilter() ,
   a_SecondaryVertexFitter(p.p_VertexPull) , 
   a_CTTest_Plane0(0) ,
   a_CTTest_Plane1(1) ,
   a_CTTest_Plane2(2) 
{
   std::cout << "Building SelectionManager" << std::endl;

   // Set the selection parameters
   TheParams = p;
   DeclareCuts();
   a_SelectorBDTManager.SetCuts(p.p_Proton_PID_Cut,p.p_Pion_PID_Cut,p.p_Separation_Cut);
   a_AnalysisBDTManager.SetPull(p.p_VertexPull);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Close selection manager

void SelectionManager::Close(){
   if(f_Hists != nullptr) f_Hists->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Set the beam being used

void SelectionManager::SetBeamMode(int mode){

   assert(mode == kFHC || mode == kRHC || mode == kBNB);

   BeamMode = mode;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Return the selection parameters loaded

SelectionParameters SelectionManager::GetParams(){
   return TheParams;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load a sample of events, specify name, type and POT. Can also load list of selected events
// if using the event list cut

void SelectionManager::AddSample(std::string Name,std::string Type,double SamplePOT,std::string EventList){

   std::cout << "Processing Sample " << Name << " of type " << Type << " and POT " << SamplePOT <<  std::endl;

   if(Type != "Data") thisSampleWeight = POT/SamplePOT;
   else thisSampleWeight = 1.0;

   thisSampleName = Name;
   thisSampleType = Type;

   if(thisSampleType != "Data" && thisSampleType != "EXT"){
      if (thisSampleName.find("GENIE") != std::string::npos) thisSampleGenerator = kGENIE;
      else if (thisSampleName.find("NuWro") != std::string::npos) thisSampleGenerator = kNuWro;      
      else thisSampleGenerator = -1;
   }

   if(Type == "Data") fHasData = true;

   TString DataDir = getenv("DATA_DIR");         

   TString file = DataDir + TString("/EventLists/") +  TString(EventList);  

   if(EventList != ""){
      std::cout << "Adding event list " << EventList << std::endl;
      a_EventListFilter.SetList(file); 
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load an event, normalise its weight using the POT of the current sample

void SelectionManager::AddEvent(Event &e){

   // Sample Orthogonality
   //if(thisSampleType == "Hyperon" && !e.EventHasHyperon){ e.Weight = 0.0; return; }
   //if(thisSampleType == "Background" && e.EventHasHyperon){ e.Weight = 0.0; return; }

   if(thisSampleType == "Neutron" && !e.EventHasNeutronScatter){ e.Weight = 0.0; return; }
   if(thisSampleType != "Neutron" && e.EventHasNeutronScatter){ e.Weight = 0.0; return; }
   if(thisSampleType != "Hyperon" &&  e.EventHasHyperon){ e.Weight = 0.0; return; }
   if(thisSampleType == "Hyperon" &&  !e.EventHasHyperon){ e.Weight = 0.0; return; }

   // Set flux weight if setup
   if(thisSampleType != "Data" && thisSampleType != "EXT"){
      if(fUseFluxWeight) e.Weight *= a_FluxWeightCalc.GetFluxWeight(e);        
      if(fUseGenWeight) {
         a_GenG4WeightCalc.LoadEvent(e);
         e.Weight *= a_GenG4WeightCalc.GetCVWeight();
      }
   }

   if(thisSampleType != "Data") e.Weight *= thisSampleWeight;

   for(size_t i_c=0;i_c<Cuts.size();i_c++){
      Cuts[i_c].fTotalEvents += e.Weight;    
      if(e.EventIsSignal) Cuts[i_c].fSignalEvents += e.Weight;
      if(e.GoodReco) Cuts[i_c].fGoodRecoEvents += e.Weight;

      Cuts[i_c].fTotalEventsVar += e.Weight*e.Weight;    
      if(e.EventIsSignal) Cuts[i_c].fSignalEventsVar += e.Weight*e.Weight;
      if(e.GoodReco) Cuts[i_c].fGoodRecoEventsVar += e.Weight*e.Weight;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void SelectionManager::UseFluxWeight(bool usefluxweight){
   fUseFluxWeight = usefluxweight;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void SelectionManager::UseGenWeight(bool usegenweight){
   fUseGenWeight = usegenweight;
} 

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the signal definition

void SelectionManager::SetSignal(Event &e){

   e.EventIsSignal = false;
   e.EventIsSignalSigmaZero = false;
   e.GoodReco = false;

   std::vector<bool> IsSignal_tmp = e.IsSignal;
   std::vector<bool> IsSignalSigmaZero_tmp = e.IsSignalSigmaZero;

   for(size_t i_tr=0;i_tr<e.NMCTruths;i_tr++){

      IsSignal_tmp.at(i_tr) = false;
      IsSignalSigmaZero_tmp.at(i_tr) = false;

      e.InActiveTPC.at(i_tr) = a_FiducialVolume.InFiducialVolume(e.TruePrimaryVertex.at(i_tr)); 

      if(e.IsSignal.at(i_tr) || e.IsSignalSigmaZero.at(i_tr)){

         bool found_proton=false,found_pion=false;

         for(size_t i_d=0;i_d<e.Decay.size();i_d++){

            if(e.Decay.at(i_d).MCTruthIndex == i_tr && e.Decay.at(i_d).PDG == 2212 && e.Decay.at(i_d).ModMomentum > 0.3) 
               found_proton = true;

            if(e.Decay.at(i_d).MCTruthIndex == i_tr && e.Decay.at(i_d).PDG == -211 && e.Decay.at(i_d).ModMomentum > 0.1) 
               found_pion = true;
         }                   

         IsSignal_tmp.at(i_tr) = found_proton && found_pion && e.InActiveTPC.at(i_tr) && e.IsSignal.at(i_tr);
         IsSignalSigmaZero_tmp.at(i_tr) = found_proton && found_pion && e.InActiveTPC.at(i_tr) && e.IsSignalSigmaZero.at(i_tr);
      }
   }

   e.IsSignal = IsSignal_tmp;
   e.IsSignalSigmaZero = IsSignalSigmaZero_tmp;

   e.EventIsSignal = std::find(e.IsSignal.begin(),e.IsSignal.end(), true) != e.IsSignal.end();
   e.EventIsSignalSigmaZero = std::find(e.IsSignalSigmaZero.begin(),e.IsSignalSigmaZero.end(), true) != e.IsSignalSigmaZero.end();

   // Search the list of reco'd tracks for the proton and pion
   bool found_proton=false,found_pion=false;

   if(!e.EventIsSignal) return;

   for(size_t i_tr=0;i_tr<e.TracklikePrimaryDaughters.size();i_tr++){
      if(e.TracklikePrimaryDaughters.at(i_tr).HasTruth && e.TracklikePrimaryDaughters.at(i_tr).TrackTruePDG == 2212 && e.TracklikePrimaryDaughters.at(i_tr).TrackTrueOrigin == 2){ 
         found_proton = true; 
         e.TrueDecayProtonIndex = i_tr; 
      }
      if(e.TracklikePrimaryDaughters.at(i_tr).HasTruth && e.TracklikePrimaryDaughters.at(i_tr).TrackTruePDG == -211 && e.TracklikePrimaryDaughters.at(i_tr).TrackTrueOrigin == 2){
         found_pion = true;
         e.TrueDecayPionIndex = i_tr; 
      }
   }

   e.GoodReco = e.EventIsSignal && found_proton && found_pion; 
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load the list of cuts

void SelectionManager::DeclareCuts(){

   for(size_t i_c=0;i_c<CutTypes.size();i_c++){
      Cut c;
      c.fName = CutTypes.at(i_c);
      Cuts.push_back(c);
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Add an event to a cut

void SelectionManager::UpdateCut(Event e,bool Passed,std::string CutName){

   for(size_t i_c=0;i_c<Cuts.size();i_c++){

      if(Cuts.at(i_c).fName == CutName){

         Cuts[i_c].fEventsIn += e.Weight;
         if(e.EventIsSignal) Cuts[i_c].fSignalEventsIn += e.Weight;
         if(e.GoodReco) Cuts[i_c].fGoodRecoEventsIn += e.Weight;

         Cuts[i_c].fEventsInVar += e.Weight*e.Weight;
         if(e.EventIsSignal) Cuts[i_c].fSignalEventsInVar += e.Weight*e.Weight;
         if(e.GoodReco) Cuts[i_c].fGoodRecoEventsInVar += e.Weight*e.Weight;

         if(Passed){

            Cuts[i_c].fEventsOut += e.Weight;
            if(e.EventIsSignal) Cuts[i_c].fSignalEventsOut += e.Weight;
            if(e.GoodReco) Cuts[i_c].fGoodRecoEventsOut += e.Weight;

            Cuts[i_c].fEventsOutVar += e.Weight*e.Weight;
            if(e.EventIsSignal) Cuts[i_c].fSignalEventsOutVar += e.Weight*e.Weight;
            if(e.GoodReco) Cuts[i_c].fGoodRecoEventsOutVar += e.Weight*e.Weight;
         }       
      }
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Retrieve the results from a cut

Cut SelectionManager::GetCut(std::string CutName){

   for(size_t i_c=0;i_c<Cuts.size();i_c++){

      if(Cuts.at(i_c).fName == CutName) return Cuts.at(i_c);

   }

   std::cout << "Cut " << CutName << " not found" << std::endl;

   Cut c;
   return c;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Reset the data in the cuts

void SelectionManager::Reset(){
   for(size_t i_c=0;i_c<Cuts.size();i_c++) Cuts[i_c].Reset();  
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load the track selector BDT weights

void SelectionManager::ImportSelectorBDTWeights(std::string WeightDir){

   std::cout << "SelectionManager: Importing Selector BDT Weights from " << WeightDir << std::endl;
   a_SelectorBDTManager.SetupSelectorBDT(WeightDir);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Load the analysis BDT weights

void SelectionManager::ImportAnalysisBDTWeights(std::string WeightDir){

   std::cout << "SelectionManager: Importing Analysis BDT Weights from " << WeightDir << std::endl;
   a_AnalysisBDTManager.SetupAnalysisBDT(WeightDir);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the fiducial volume cut

bool SelectionManager::FiducialVolumeCut(Event e){
   bool passed = a_FiducialVolume.InFiducialVolume(e.RecoPrimaryVertex);

   UpdateCut(e,passed,"FV");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the three track cut

bool SelectionManager::TrackCut(Event e){

   bool passed = e.NPrimaryTrackDaughters > 2; 

   UpdateCut(e,passed,"Tracks");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the zero shower cut

bool SelectionManager::ShowerCut(Event e){

   bool passed = e.NPrimaryShowerDaughters < 1; 

   UpdateCut(e,passed,"Showers");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Select the muon candidate

bool SelectionManager::ChooseMuonCandidate(Event &e){

   int i_muon = a_MuonID.SelectCandidate(e.TracklikePrimaryDaughters);

   if(i_muon == -1){
      UpdateCut(e,false,"MuonID");
      return false;
   }

   RecoParticle theMuonCandidate = e.TracklikePrimaryDaughters.at(i_muon);
   e.MuonCandidate = theMuonCandidate;
   e.TracklikePrimaryDaughters.erase(e.TracklikePrimaryDaughters.begin()+i_muon);

   UpdateCut(e,true,"MuonID");

   return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the secondary track length cut

bool SelectionManager::TrackLengthCut(Event e){

   bool passed = a_TrackLengthCutManager.ApplyCut(e.TracklikePrimaryDaughters);

   UpdateCut(e,passed,"SubleadingTracks");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Select proton and pion tracks. If cheat is selected uses truth information instead of BDT 

bool SelectionManager::ChooseProtonPionCandidates(Event &e,bool cheat){

   std::pair<int,int> candidates;

   if(cheat) candidates = a_SelectorBDTManager.NominateTracksCheat(e);
   else candidates = a_SelectorBDTManager.NominateTracks(e);

   bool passed = candidates.first != -1 && candidates.second != -1;

   if(!passed){
      UpdateCut(e,passed,"DecaySelector");
      return false;
   }

   e.DecayProtonCandidate = e.TracklikePrimaryDaughters.at(candidates.first);
   e.DecayPionCandidate = e.TracklikePrimaryDaughters.at(candidates.second);

   // Erase from track vector
   std::vector<RecoParticle> TracklikePrimaryDaughters_tmp;

   for(size_t i=0;i<e.TracklikePrimaryDaughters.size();i++)
      if(i != candidates.first && i != candidates.second) TracklikePrimaryDaughters_tmp.push_back(e.TracklikePrimaryDaughters.at(i));

   e.TracklikePrimaryDaughters = TracklikePrimaryDaughters_tmp;

   UpdateCut(e,passed,"DecaySelector");
   return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Calculate decay analysis BDT score and apply cut

bool SelectionManager::AnalysisBDTCut(Event &e){

   bool passed = a_AnalysisBDTManager.CalculateScore(e) > TheParams.p_AnalysisBDT_Cut; 

   UpdateCut(e,passed,"DecayAnalysis");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Select events in the list already loaded

bool SelectionManager::EventListCut(Event e){

   bool passed = a_EventListFilter.EventPassed(e.run,e.subrun,e.event); 

   UpdateCut(e,passed,"Connectedness");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the connectedness test

bool SelectionManager::ConnectednessTest(Event e, int nplanes){

   int muon_index = e.MuonCandidate.Index;
   int proton_index = e.DecayProtonCandidate.Index;
   int pion_index = e.DecayPionCandidate.Index;

   a_CTTest_Plane0.LoadInfo(e.ConnSeedIndexes_Plane0,e.ConnOutputIndexes_Plane0,e.ConnOutputSizes_Plane0,e.ConnSeedChannels_Plane0);
   a_CTTest_Plane1.LoadInfo(e.ConnSeedIndexes_Plane1,e.ConnOutputIndexes_Plane1,e.ConnOutputSizes_Plane1,e.ConnSeedChannels_Plane1);
   a_CTTest_Plane2.LoadInfo(e.ConnSeedIndexes_Plane2,e.ConnOutputIndexes_Plane2,e.ConnOutputSizes_Plane2,e.ConnSeedChannels_Plane2);

   int npassed = 0;
   if(a_CTTest_Plane0.DoTest(muon_index,proton_index,pion_index)) npassed++;
   if(a_CTTest_Plane1.DoTest(muon_index,proton_index,pion_index)) npassed++;
   if(a_CTTest_Plane2.DoTest(muon_index,proton_index,pion_index)) npassed++;


   bool passed = npassed >= nplanes; 

   UpdateCut(e,passed,"Connectedness");

   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the invariant mass cut

bool SelectionManager::WCut(Event e){

   double W = ProtonPionInvariantMass(e.DecayProtonCandidate,e.DecayPionCandidate); 
   bool passed = W > TheParams.p_W_Min && W < TheParams.p_W_Max;
   UpdateCut(e,passed,"InvariantMass");
   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Apply the angle cut

bool SelectionManager::AngleCut(Event e){

   SecondaryVertex V = a_SecondaryVertexFitter.MakeVertex(e.DecayProtonCandidate,e.DecayPionCandidate);
   TVector3 GapVector = V.Vertex - e.RecoPrimaryVertex;
   double Reco_DecayLength = (V.Vertex - e.RecoPrimaryVertex).Mag();
   TLorentzVector Lambda4Momentum = ProtonPion4Momentum(e.DecayProtonCandidate,e.DecayPionCandidate);
   TVector3 LambdaDirection(Lambda4Momentum.X(),Lambda4Momentum.Y(),Lambda4Momentum.Z());
   double alpha = (180/3.1415)*LambdaDirection.Angle(GapVector);

   bool passed = alpha < TheParams.p_Alpha_Cut;
   UpdateCut(e,passed,"AlphaAngle");
   return passed;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Setup histograms, vector with bin boundaries and string containing axis titles (same format
// as Root TH1 titles).

void SelectionManager::SetupHistograms(std::vector<double> boundaries,std::string title){

   if(boundaries.size() < 2)
      throw std::invalid_argument("Require at least 2 bin boundaries"); 

   fTitle = title;
   fHistNBins = boundaries.size()-1;
   fHistLow = boundaries.front();
   fHistHigh = boundaries.back();
   fHistBoundaries = boundaries;

   const int arr_n = boundaries.size();
   Double_t arr_boundaries[arr_n];
   for(size_t i=0;i<arr_n;i++) arr_boundaries[i] = boundaries.at(i);

   for(size_t i_proc=0;i_proc<Procs.size();i_proc++){
      std::string histname = "h_" + Procs.at(i_proc);
      Hists_ByProc[Procs.at(i_proc)] = new TH1D(histname.c_str(),fTitle.c_str(),fHistNBins,arr_boundaries);
   }

   for(size_t i_type=0;i_type<Types.size();i_type++){
      std::string histname = "h_ByType_" + Types.at(i_type);
      Hists_ByType[Types.at(i_type)] = new TH1D(histname.c_str(),fTitle.c_str(),fHistNBins,arr_boundaries);
   }

   for(size_t i_type=0;i_type<Types2.size();i_type++){
      std::string histname = "h_ByType2_" + Types2.at(i_type);
      Hists_ByType2[Types2.at(i_type)] = new TH1D(histname.c_str(),fTitle.c_str(),fHistNBins,arr_boundaries);
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Setup histograms. Give bin boundaries and number of bins. String containing axis titles 
// (same format as Root TH1 titles).

void SelectionManager::SetupHistograms(int n,double low,double high,std::string title){

   fTitle = title;
   fHistNBins = n;
   fHistLow = low;
   fHistHigh = high;

   for(size_t i_proc=0;i_proc<Procs.size();i_proc++){
      std::string histname = "h_" + Procs.at(i_proc);
      Hists_ByProc[Procs.at(i_proc)] = new TH1D(histname.c_str(),fTitle.c_str(),n,low,high);
   }

   for(size_t i_type=0;i_type<Types.size();i_type++){
      std::string histname = "h_ByType_" + Types.at(i_type);
      Hists_ByType[Types.at(i_type)] = new TH1D(histname.c_str(),fTitle.c_str(),n,low,high);
   }

   for(size_t i_type=0;i_type<Types2.size();i_type++){
      std::string histname = "h_ByType2_" + Types2.at(i_type);
      Hists_ByType2[Types2.at(i_type)] = new TH1D(histname.c_str(),fTitle.c_str(),n,low,high);
   }

   double width = (high-low)/n;
   for(int i=0;i<n+1;i++) fHistBoundaries.push_back(low+width*i); 
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Setup histograms for a systematic. Give the type (multisim, single unisim or dual unisim),
// the number of universes and name

void SelectionManager::AddSystematic(int systype,int universes,std::string name){

   std::cout << "Setting up systematic " << name << std::endl;

   const int arr_n = fHistBoundaries.size();
   Double_t arr_boundaries[arr_n];
   for(size_t i=0;i<arr_n;i++) arr_boundaries[i] = fHistBoundaries.at(i);

   if(systype == kMultisim){
      for(size_t i_type=0;i_type<Types.size();i_type++){
         Multisim_Sys_Hists[Types.at(i_type)][name].resize(universes);
         for(size_t i=0;i<universes;i++) 
            Multisim_Sys_Hists[Types.at(i_type)][name].at(i) = new TH1D(("h_" + Types.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);      
      }
      for(size_t i_type=0;i_type<Types2.size()-4;i_type++){
         Multisim_Sys_Hists[Types2.at(i_type)][name].resize(universes);
         for(size_t i=0;i<universes;i++) 
            Multisim_Sys_Hists[Types2.at(i_type)][name].at(i) = new TH1D(("h_" + Types2.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);   
      }
   }
   else if(systype == kSingleUnisim){
      for(size_t i_type=0;i_type<Types.size();i_type++){
         SingleUnisim_Sys_Hists[Types.at(i_type)][name].resize(1);
         for(size_t i=0;i<1;i++)
            SingleUnisim_Sys_Hists[Types.at(i_type)][name][i] = new TH1D(("h_" + Types.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);    
      }
      for(size_t i_type=0;i_type<Types2.size()-4;i_type++){
         SingleUnisim_Sys_Hists[Types2.at(i_type)][name].resize(1);
         for(size_t i=0;i<1;i++)
            SingleUnisim_Sys_Hists[Types2.at(i_type)][name][i] = new TH1D(("h_" + Types2.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);  
      }
   }
   else if(systype == kDualUnisim){
      for(size_t i_type=0;i_type<Types.size();i_type++){
         DualUnisim_Sys_Hists[Types.at(i_type)][name].resize(2);
         for(size_t i=0;i<2;i++) 
            DualUnisim_Sys_Hists[Types.at(i_type)][name].at(i) = new TH1D(("h_" + Types.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);
      }
      for(size_t i_type=0;i_type<Types2.size()-4;i_type++){
         DualUnisim_Sys_Hists[Types2.at(i_type)][name].resize(2);
         for(size_t i=0;i<2;i++) 
            DualUnisim_Sys_Hists[Types2.at(i_type)][name].at(i) = new TH1D(("h_" + Types2.at(i_type) + "_" + name + "_u_" + std::to_string(i)).c_str(),"",fHistNBins,arr_boundaries);
      }
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the event category 

std::string SelectionManager::GetMode(Event e){

   std::string mode;
   if(thisSampleType == "Data") mode = "Data";
   else if(thisSampleType == "EXT") mode = "EXT";
   else if(thisSampleType == "Dirt") mode = "Dirt";
   else if(e.EventIsSignal) mode = "Signal";
   else if(e.Mode.at(0) == "HYP") mode = "OtherHYP";
   else  mode = "OtherNu";

   return mode;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the event category using different labelling convention

std::string SelectionManager::GetMode2(Event e){

   std::string mode;

   if(thisSampleType == "Data") return "Data";
   else if(thisSampleType == "EXT") return "EXT";
   else if (thisSampleType == "Dirt") return "Dirt";

   bool islambdacharged =  std::find(e.IsLambdaCharged.begin(),e.IsLambdaCharged.end(), true) != e.IsLambdaCharged.end();
   bool ishyperon = e.Hyperon.size(); 

   if(e.EventIsSignal) return "DirectLambda";
   else if(e.Mode.at(0) == "HYP") return "DirectHYP"; 
   else if(e.Mode.at(0) == "RES" && islambdacharged) return "RESLambda";
   else if(e.Mode.at(0) == "RES" && e.Hyperon.size()) return "RESHYP"; 
   else if(e.Mode.at(0) == "DIS" && islambdacharged) return "DISLambda";
   else if(e.Mode.at(0) == "DIS" && e.Hyperon.size()) return "DISHYP"; 
   else if(e.EventHasNeutronScatter) return "Neutron";
   else return "Other";

   return mode;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Fill the histograms

void SelectionManager::FillHistograms(Event e,double variable,double weight){

   std::string mode;

   bool isNuBackground = false;

   if( thisSampleType == "Data" ) mode = "Data";
   else if( thisSampleType == "EXT" ) mode = "EXT";
   else if( thisSampleType == "Dirt" ) mode = "Dirt";
   else  if( e.EventIsSignal ) mode = "Signal";
   else if( e.Mode.at(0) == "HYP") mode = "OtherHYP";
   else { mode = e.Mode.at(0); isNuBackground = true; }

   if(!isNuBackground){
      Hists_ByType[ mode ]->Fill(variable,weight*e.Weight);
      Hists_ByProc[ mode ]->Fill(variable,weight*e.Weight);
   }
   else {

      Hists_ByType[ "OtherNu" ]->Fill(variable,weight*e.Weight);

      if(mode != "ElectronScattering" && mode != "Diffractive" && mode != "Other"){
         if(e.CCNC.at(0) == "CC"){
            Hists_ByProc[ "CC"+mode ]->Fill(variable,weight*e.Weight);
         }
         else Hists_ByProc[ "NC" ]->Fill(variable,weight*e.Weight);
      }
      else if( mode == "ElectronScattering") Hists_ByProc[ "ElectronScattering" ]->Fill(variable,weight*e.Weight);
      else if( mode == "Diffractive") Hists_ByProc[ "Diffractive" ]->Fill(variable,weight*e.Weight);
      else Hists_ByProc[ "Other" ]->Fill(variable,weight*e.Weight);

   }

   if(mode != "Data"){
      Hists_ByType["All"]->Fill(variable,weight*e.Weight);
      Hists_ByProc["All"]->Fill(variable,weight*e.Weight);
   }


   std::string mode2 = GetMode2(e);
   Hists_ByType2[mode2]->Fill(variable,weight*e.Weight);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Fill systematics histogram for a single universe with a single weight

void SelectionManager::FillHistogramsSys(Event e,double variable,std::string name,int universe,double weight){

   if(thisSampleType == "Data") return;

   if(std::isnan(e.Weight) || std::isnan(weight) || std::isinf(e.Weight) || std::isinf(weight)){
       //std::cout << "Nan weight detected for dial " << name << "  event " << e.run << " " << e.subrun << " " << e.event << " skipping" << std::endl;
      return;
   }

   std::string mode = GetMode(e);
   std::string mode2 = GetMode2(e);

   if (Multisim_Sys_Hists[mode].find(name) != Multisim_Sys_Hists[mode].end()){
      Multisim_Sys_Hists[mode][name].at(universe)->Fill(variable,weight*e.Weight);
      Multisim_Sys_Hists["All"][name].at(universe)->Fill(variable,weight*e.Weight);
      if(mode2 != "Dirt" && mode2 != "EXT" && mode2 != "Data") Multisim_Sys_Hists[mode2][name].at(universe)->Fill(variable,weight*e.Weight);
      return;
   }

   else if (SingleUnisim_Sys_Hists[mode].find(name) != SingleUnisim_Sys_Hists[mode].end()){
      SingleUnisim_Sys_Hists[mode][name].at(universe)->Fill(variable,weight*e.Weight);
      SingleUnisim_Sys_Hists["All"][name].at(universe)->Fill(variable,weight*e.Weight);
      if(mode2 != "Dirt" && mode2 != "EXT" && mode2 != "Data") SingleUnisim_Sys_Hists[mode2][name].at(universe)->Fill(variable,weight*e.Weight);
      return;
   }

   else if (DualUnisim_Sys_Hists[mode].find(name) != DualUnisim_Sys_Hists[mode].end()){
      DualUnisim_Sys_Hists[mode][name].at(universe)->Fill(variable,weight*e.Weight);
      DualUnisim_Sys_Hists["All"][name].at(universe)->Fill(variable,weight*e.Weight);
      if(mode2 != "Dirt" && mode2 != "EXT" && mode2 != "Data") DualUnisim_Sys_Hists[mode2][name].at(universe)->Fill(variable,weight*e.Weight);
      return;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Fill systematics histogram for all universes. Supply vector containing weights for different
// universes

void SelectionManager::FillHistogramsSys(Event e,double variable,std::string name,std::vector<double> weights){

   for(size_t i=0;i<weights.size();i++)
      FillHistogramsSys(e,variable,name,i,weights.at(i)); 

}

///////////////////////////////////////////////////////////////////////////////////////////////

// Draw the histograms

void SelectionManager::DrawHistograms(std::string label,double Scale,double SignalScale){

   OpenHistFile(label);

   double y_limit = -1;
   if(y_limit != -1) gStyle->SetHistTopMargin(0);

   system(("mkdir -p " + PlotDir).c_str());

   // Create weight sums

   for(size_t i_proc=0;i_proc<Procs.size();i_proc++){
      Hists_ByProc[Procs.at(i_proc)]->Scale(Scale);
      Hists_ByProc[Procs.at(i_proc)]->Sumw2();
   }

   for(size_t i_type=0;i_type<Types.size();i_type++){
      Hists_ByType[Types.at(i_type)]->Scale(Scale);
      Hists_ByType[Types.at(i_type)]->Sumw2();
   }

   for(size_t i_type=0;i_type<Types2.size();i_type++){
      Hists_ByType2[Types2.at(i_type)]->Scale(Scale);
      Hists_ByType2[Types2.at(i_type)]->Sumw2();
   }


   Hists_ByProc["Signal"]->Scale(SignalScale);
   Hists_ByType["Signal"]->Scale(SignalScale);

   TH1D* h_errors = (TH1D*)Hists_ByType["All"]->Clone("h_errors");

   std::vector<std::string> captions = {"Signal","Other HYP","Other #nu","Dirt","EXT","Data"};
   std::vector<int> colors = {8,46,38,30,15,0};
   std::vector<TH1D*> Hists_ByType_v = {Hists_ByType["Signal"],Hists_ByType["OtherHYP"],Hists_ByType["OtherNu"],Hists_ByType["Dirt"],Hists_ByType["EXT"],Hists_ByType["Data"]};
   DrawHistogram(Hists_ByType_v,h_errors,captions,PlotDir,label+"_ByType",BeamMode,Run,POT,SignalScale,fHasData,colors,BinLabels,std::make_pair(0,0));

   std::vector<std::string> captions2 = {"Direct #Lambda","Direct Hyp","Neutron","Dirt","RES #Lambda","RES Hyp","Other #nu","EXT","DIS #Lambda","DIS Hyp","Data"};
   std::vector<int> colors2 = {8,kBlue-7,kRed-7,kCyan+3,kGreen+3,kBlue-10,kRed-10,kMagenta-7,30,15,0};

   std::vector<TH1D*> Hists_ByType_v2 = {Hists_ByType2["DirectLambda"],Hists_ByType2["DirectHYP"],Hists_ByType2["Neutron"],Hists_ByType2["Dirt"],Hists_ByType2["RESLambda"],Hists_ByType2["RESHYP"],Hists_ByType2["Other"],Hists_ByType2["EXT"],Hists_ByType2["DISLambda"],Hists_ByType2["DISHYP"],Hists_ByType2["Data"]};

   DrawHistogram(Hists_ByType_v2,h_errors,captions2,PlotDir,label+"_ByType2",BeamMode,Run,POT,SignalScale,fHasData,colors2,BinLabels,std::make_pair(0,0));

   std::vector<std::string> captions3 = {"Signal","Other HYP","CCQEL","CCRES","CCDIS","CCMEC","CCCOH","NC","ElectronScattering","Diffractive","Other","Dirt","EXT","Data"};
   std::vector<int> colors3 = {8,46,2,3,4,5,6,7,9,11,12,30,15,0};
   std::vector<TH1D*> Hists_ByProc_v = {Hists_ByProc["Signal"],Hists_ByProc["OtherHYP"],Hists_ByProc["CCQEL"],Hists_ByProc["CCRES"],Hists_ByProc["CCDIS"],Hists_ByProc["CCMEC"],Hists_ByProc["CCCOH"],Hists_ByProc["NC"],Hists_ByProc["ElectronScattering"],Hists_ByProc["Diffractive"],Hists_ByProc["Other"],Hists_ByProc["EXT"],Hists_ByProc["Dirt"],Hists_ByProc["Data"]};
   DrawHistogram(Hists_ByProc_v,h_errors,captions3,PlotDir,label+"_ByProc",BeamMode,Run,POT,SignalScale,fHasData,colors3,BinLabels,std::make_pair(0,0));

   std::map<std::string,TH1D*>::iterator it;
   for (it = Hists_ByType.begin(); it != Hists_ByType.end(); it++)
      it->second->Write(it->first.c_str());
   for (it = Hists_ByType2.begin(); it != Hists_ByType2.end(); it++){
      if(it->first == "Data" || it->first == "All" || it->first == "EXT" || it->first == "Dirt") continue;
      it->second->Write(it->first.c_str());
   }

   h_errors->Write("ErrorBand");
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Draw the systematics histograms

void SelectionManager::DrawHistogramsSys(std::string label,std::string name,std::string type){

   std::cout << "Drawing sys histograms label=" << label << "  name=" << name << "  type=" << type << std::endl;

   system(("mkdir -p " + PlotDir).c_str());
   system(("mkdir -p " + RootfileDir).c_str());

   int systype=-1;
   std::vector<TH1D*> ToDraw;

   if (Multisim_Sys_Hists[type].find(name) != Multisim_Sys_Hists[type].end()){
      ToDraw = Multisim_Sys_Hists[type][name]; 
      systype = kMultisim;
   }
   else if (SingleUnisim_Sys_Hists[type].find(name) != SingleUnisim_Sys_Hists[type].end()){
      ToDraw = SingleUnisim_Sys_Hists[type][name];
      systype = kSingleUnisim;
   }
   else if (DualUnisim_Sys_Hists[type].find(name) != DualUnisim_Sys_Hists[type].end()){
      ToDraw = DualUnisim_Sys_Hists[type][name];
      systype = kDualUnisim;
   }

   if(!ToDraw.size()){
      std::cout << "Can't find systematic name=" << name << " type=" << type << std::endl;
      return;
   }

   if(BinLabels.size())
      for(size_t i_u=0;i_u<ToDraw.size();i_u++)
         for(int i=1;i<fHistNBins+1;i++) ToDraw.at(i_u)->GetXaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());

   THStack *hs = new THStack("hs",fTitle.c_str());

   TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
   l->SetBorderSize(0);
   l->SetNColumns(2);
   l->AddEntry(Hists_ByType[type],"Central Value","L");

   // Get the maximum of all the histograms
   double maximum = 0.0;

   f_Hists->mkdir(name.c_str()); 
   f_Hists->GetDirectory(name.c_str())->cd();

   for(size_t i=0;i<ToDraw.size();i++){

      if(systype == kMultisim){
         ToDraw.at(i)->SetLineColor(kGreen);
         if(i == 0) l->AddEntry(ToDraw.at(0),"Variations","L");
      }
      else if(systype == kSingleUnisim){
         ToDraw.at(i)->SetLineColor(kRed);
         l->AddEntry(ToDraw.at(0),"Alternative Model","L");
      }
      else if(systype == kDualUnisim && i == 0){
         ToDraw.at(i)->SetLineColor(kRed);
         l->AddEntry(ToDraw.at(i),"+ 1 #sigma","L");
      }
      else if(systype == kDualUnisim && i == 1){
         ToDraw.at(i)->SetLineColor(kBlue); 
         l->AddEntry(ToDraw.at(i),"- 1 #sigma","L");
      }

      hs->Add(ToDraw.at(i));

      if(GetHistMax(ToDraw.at(i)) > maximum) maximum = GetHistMax(ToDraw.at(i));

      ToDraw.at(i)->SetDirectory(new TDirectory(name.c_str(),name.c_str()));      
      ToDraw.at(i)->Write();
   }

   f_Hists->cd();

   TH1D* h_CV = nullptr;
   if(std::find(Types.begin(), Types.end(), type) != Types.end())
      h_CV = Hists_ByType[type];
   else if(std::find(Types2.begin(), Types2.end(), type) != Types2.end())
      h_CV = Hists_ByType2[type];   
   else {
      std::cout << "Can't find histogram for channel=" << type << std::endl;
      return;
   }

   if(GetHistMax(h_CV) > maximum) maximum = GetHistMax(h_CV);
   h_CV->SetLineColor(1);
   h_CV->SetLineWidth(2);
   h_CV->SetFillColor(0);
   h_CV->SetMaximum(maximum*1.25);
   hs->Add(h_CV);

   TCanvas *c = new TCanvas("c","c",800,600);

   TPad *p_plot = new TPad("pad1","pad1",0,0,1,0.85);
   TPad *p_legend = new TPad("pad2","pad2",0,0.85,1,1);

   p_legend->SetBottomMargin(0);
   p_legend->SetTopMargin(0.1);
   p_plot->SetTopMargin(0.01);

   TLegend *l_Watermark = new TLegend(0.45,0.900,0.89,0.985);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.005);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetTextFont(62);
   l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");
   //l_Watermark->SetHeader("MicroBooNE Simulation In Progress","R");

   TLegend *l_POT = new TLegend(0.55,0.820,0.89,0.900);
   l_POT->SetBorderSize(0);
   l_POT->SetMargin(0.005);
   l_POT->SetTextAlign(32);
   l_POT->SetTextSize(0.05);
   if(BeamMode == kFHC) l_POT->SetHeader(("NuMI FHC, " + to_string_with_precision(POT/1e20,1) + " #times 10^{20} POT").c_str());
   if(BeamMode == kRHC) l_POT->SetHeader(("NuMI RHC, " + to_string_with_precision(POT/1e20,1) + " #times 10^{20} POT").c_str());
   if(BeamMode == kBNB) l_POT->SetHeader(("BNB, " + to_string_with_precision(POT/1e20,1) + " #times 10^{20} POT").c_str());

   p_legend->Draw();
   p_legend->cd();
   l->Draw();
   c->cd();
   p_plot->Draw();
   p_plot->cd();

   hs->Draw("HIST nostack");
   hs->GetXaxis()->SetTitleSize(0.05);
   hs->GetYaxis()->SetTitleSize(0.05);
   hs->GetXaxis()->SetTitleOffset(0.9);
   hs->GetYaxis()->SetTitleOffset(0.9);

   hs->GetXaxis()->SetLabelSize(0.045);
   hs->GetYaxis()->SetLabelSize(0.045);

   h_CV->Draw("HIST same");

   if(BinLabels.size()){
      for(int i=1;i<fHistNBins+1;i++) hs->GetXaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());
      hs->GetXaxis()->SetLabelSize(0.07);
   }

   p_plot->Update(); 

   if(POT > 0) l_POT->Draw();
   l_Watermark->Draw();

   c->Print((PlotDir + label + "_Sys_" + type + "_" +  name + ".png").c_str());
   c->Print((PlotDir + label + "_Sys_" + type + "_" +  name + ".pdf").c_str());
   c->Print((PlotDir + label + "_Sys_" + type + "_" +  name + ".C").c_str());

   c->Close();

   std::cout << "Done drawing systematic hists for dial " << name << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the covariace matrix for systematic dial

TMatrixD SelectionManager::GetCovarianceMatrix(std::string label,std::string name,std::string type){

   system(("mkdir -p " + PlotDir).c_str());

   int systype=-1;

   std::cout << "Getting Covariance Matrix for " << name << std::endl;

   OpenHistFile(label);

   std::vector<TH1D*> ToUse;

   if (Multisim_Sys_Hists[type].find(name) != Multisim_Sys_Hists[type].end()){
      ToUse = Multisim_Sys_Hists[type][name];
      systype = kMultisim;
   }
   else if (SingleUnisim_Sys_Hists[type].find(name) != SingleUnisim_Sys_Hists[type].end()){
      ToUse = SingleUnisim_Sys_Hists[type][name];
      systype = kSingleUnisim;
   }
   else if (DualUnisim_Sys_Hists[type].find(name) != DualUnisim_Sys_Hists[type].end()){
      ToUse = DualUnisim_Sys_Hists[type][name];
      systype = kDualUnisim;
   }

   const int arr_n = fHistBoundaries.size();
   Double_t arr_boundaries[arr_n];
   for(size_t i=0;i<arr_n;i++) arr_boundaries[i] = fHistBoundaries.at(i);

   std::string title = ";" + std::string(Hists_ByType["All"]->GetXaxis()->GetTitle()) + ";" + std::string(Hists_ByType["All"]->GetXaxis()->GetTitle());

   TMatrixD Cov(fHistNBins,fHistNBins);
   TH2D *h_Cov = new TH2D("h_Cov",title.c_str(),fHistNBins,arr_boundaries,fHistNBins,arr_boundaries);

   TMatrixD frac_Cov(fHistNBins,fHistNBins);
   TH2D *h_frac_Cov = new TH2D("h_frac_Cov",title.c_str(),fHistNBins,arr_boundaries,fHistNBins,arr_boundaries);

   if(BinLabels.size()){
      for(int i=1;i<fHistNBins+1;i++){
         h_Cov->GetXaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());
         h_Cov->GetYaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());
         h_frac_Cov->GetXaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());
         h_frac_Cov->GetYaxis()->SetBinLabel(i,BinLabels.at(i-1).c_str());
      }

      h_Cov->GetXaxis()->SetLabelSize(0.07);
      h_Cov->GetYaxis()->SetLabelSize(0.07);
      h_frac_Cov->GetXaxis()->SetLabelSize(0.07);
      h_frac_Cov->GetYaxis()->SetLabelSize(0.07);

   }

   for(int i_b=1;i_b<fHistNBins+1;i_b++){
      for(int j_b=1;j_b<fHistNBins+1;j_b++){

         if(!(Hists_ByType[type]->GetBinContent(i_b) > 0 && Hists_ByType[type]->GetBinContent(j_b) > 0)) continue;

         Cov[i_b-1][j_b-1] = 0.0;        
         frac_Cov[i_b-1][j_b-1] = 0.0; 

         double Cov_ij = 0.0;
         double Mean_i = 0.0;
         double Mean_j = 0.0;

         if(systype == kMultisim){
            for(size_t i_u=0;i_u<ToUse.size();i_u++) Mean_i += ToUse.at(i_u)->GetBinContent(i_b)/ToUse.size();
            for(size_t i_u=0;i_u<ToUse.size();i_u++) Mean_j += ToUse.at(i_u)->GetBinContent(j_b)/ToUse.size();
            for(size_t i_u=0;i_u<ToUse.size();i_u++) Cov_ij += (ToUse.at(i_u)->GetBinContent(i_b) - Mean_i)*(ToUse.at(i_u)->GetBinContent(j_b) - Mean_j)/(ToUse.size()-1);
         }
         if(systype == kSingleUnisim){
            Mean_i = Hists_ByType[type]->GetBinContent(i_b);
            Mean_j = Hists_ByType[type]->GetBinContent(j_b);
            Cov_ij = (ToUse.at(0)->GetBinContent(i_b) - Mean_i)*(ToUse.at(0)->GetBinContent(j_b) - Mean_j);
         } 
         if(systype == kDualUnisim){
            Mean_i = Hists_ByType[type]->GetBinContent(i_b);
            Mean_j = Hists_ByType[type]->GetBinContent(j_b);
            Cov_ij = (ToUse.at(1)->GetBinContent(i_b) - ToUse.at(0)->GetBinContent(i_b))*(ToUse.at(1)->GetBinContent(j_b) - ToUse.at(0)->GetBinContent(j_b))/4;
         } 
         frac_Cov[i_b-1][j_b-1] = Cov_ij/Mean_i/Mean_j; 
         Cov[i_b-1][j_b-1] = Cov_ij;
         h_Cov->SetBinContent(i_b,j_b,Cov_ij);
         h_frac_Cov->SetBinContent(i_b,j_b,Cov_ij/Mean_i/Mean_j);
      }
   }

   std::string plottitle = PlotDir + label + "_CovMatrix_Sys_" + type + "_" + name;
   DrawMatrix(Cov,h_Cov,plottitle,BinLabels.size(),fUseText);

   plottitle = PlotDir + label + "_FCovMatrix_Sys_" + type + "_" + name;
   DrawMatrix(frac_Cov,h_frac_Cov,plottitle,BinLabels.size(),fUseText);

   h_Cov->Write(("Cov_" + type + "_" + name).c_str());
   h_frac_Cov->Write(("FCov_" + type + "_" + name).c_str());

   return Cov;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the predicted number of events in a bin for a given event category

double SelectionManager::GetPrediction(int bin,std::string type){

   if(type=="") return Hists_ByType["All"]->GetBinContent(bin);
   else if(Hists_ByType.find(type) != Hists_ByType.end()) return Hists_ByType[type]->GetBinContent(bin);
   else if(Hists_ByProc.find(type) != Hists_ByProc.end()) return Hists_ByProc[type]->GetBinContent(bin);

   std::cout << "Type/Proc " << type << " not found, returning -1" << std::endl;
   return -1;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Reopen the file containing the histograms

void SelectionManager::OpenHistFile(std::string label){

   if(f_Hists == nullptr){
      system(("mkdir -p " +  RootfileDir).c_str());
      f_Hists = TFile::Open((RootfileDir + label + "_Histograms.root").c_str(),"RECREATE");
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Set the labels on the bins

void SelectionManager::SetBinLabels(std::vector<std::string> binlabels){

   if(binlabels.size() != fHistNBins) 
      throw std::invalid_argument("Set of bin labels does not match number of bins!"); 

   BinLabels = binlabels;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Divide bins by bin widths

void SelectionManager::WidthScaleHistograms(){

   std::map<std::string,TH1D*>::iterator it1;

   for (it1=Hists_ByType.begin();it1!=Hists_ByType.end();it1++){
      for(size_t i_b=1;i_b<it1->second->GetNbinsX()+1;i_b++){        
         it1->second->SetBinContent(i_b,it1->second->GetBinContent(i_b)/it1->second->GetBinWidth(i_b));
         it1->second->SetBinError(i_b,it1->second->GetBinError(i_b)/it1->second->GetBinWidth(i_b));
      }
   }

   for (it1=Hists_ByType2.begin();it1!=Hists_ByType2.end();it1++){
      for(size_t i_b=1;i_b<it1->second->GetNbinsX()+1;i_b++){        
         it1->second->SetBinContent(i_b,it1->second->GetBinContent(i_b)/it1->second->GetBinWidth(i_b));
         it1->second->SetBinError(i_b,it1->second->GetBinError(i_b)/it1->second->GetBinWidth(i_b));
      }
   }

   for (it1=Hists_ByProc.begin();it1!=Hists_ByProc.end();it1++){
      for(size_t i_b=1;i_b<it1->second->GetNbinsX()+1;i_b++){
         it1->second->SetBinContent(i_b,it1->second->GetBinContent(i_b)/it1->second->GetBinWidth(i_b));
         it1->second->SetBinError(i_b,it1->second->GetBinError(i_b)/it1->second->GetBinWidth(i_b));
      }
   }

   std::map<std::string,std::map<std::string,std::vector<TH1D*>>>::iterator it2;

   for (it2=Multisim_Sys_Hists.begin();it2!=Multisim_Sys_Hists.end();it2++){
      std::map<std::string,std::vector<TH1D*>>::iterator it22;
      for (it22=it2->second.begin();it22!=it2->second.end();it22++){
         for(size_t i=0;i<it22->second.size();i++){
            for(size_t i_b=1;i_b<it22->second.at(i)->GetNbinsX()+1;i_b++){        
               it22->second.at(i)->SetBinContent(i_b,it22->second.at(i)->GetBinContent(i_b)/it22->second.at(i)->GetBinWidth(i_b));
               it22->second.at(i)->SetBinError(i_b,it22->second.at(i)->GetBinError(i_b)/it22->second.at(i)->GetBinWidth(i_b));
            }
         }
      }
   }


   for (it2=SingleUnisim_Sys_Hists.begin();it2!=SingleUnisim_Sys_Hists.end();it2++){
      std::map<std::string,std::vector<TH1D*>>::iterator it22;
      for (it22=it2->second.begin();it22!=it2->second.end();it22++){
         for(size_t i=0;i<it22->second.size();i++){
            for(size_t i_b=1;i_b<it22->second.at(i)->GetNbinsX()+1;i_b++){        
               it22->second.at(i)->SetBinContent(i_b,it22->second.at(i)->GetBinContent(i_b)/it22->second.at(i)->GetBinWidth(i_b));
               it22->second.at(i)->SetBinError(i_b,it22->second.at(i)->GetBinError(i_b)/it22->second.at(i)->GetBinWidth(i_b));
            }
         }
      }
   }

   for (it2=DualUnisim_Sys_Hists.begin();it2!=DualUnisim_Sys_Hists.end();it2++){
      std::map<std::string,std::vector<TH1D*>>::iterator it22;
      for (it22=it2->second.begin();it22!=it2->second.end();it22++){
         for(size_t i=0;i<it22->second.size();i++){
            for(size_t i_b=1;i_b<it22->second.at(i)->GetNbinsX()+1;i_b++){        
               it22->second.at(i)->SetBinContent(i_b,it22->second.at(i)->GetBinContent(i_b)/it22->second.at(i)->GetBinWidth(i_b));
               it22->second.at(i)->SetBinError(i_b,it22->second.at(i)->GetBinError(i_b)/it22->second.at(i)->GetBinWidth(i_b));
            }
         }
      }
   }

}

///////////////////////////////////////////////////////////////////////////////////////////////

// Retrieve the vector containing bin boundaries

std::vector<double> SelectionManager::GetBinBoundaries(){
   return fHistBoundaries; 
}

///////////////////////////////////////////////////////////////////////////////////////////////

#endif
