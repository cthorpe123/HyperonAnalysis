R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so)
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so)

#include "SelectionManager.h"
#include "EventAssembler.h"
#include "Cut.h"
#include "SelectionParameters.h"

#include "Parameters.h"
#include "SampleSets_Example.h"

   void GenerateBackgroundPlots(){

      std::string label = "test";
      double POT = 1.0e21; // POT to scale samples to

      BuildTunes();

      SampleNames.push_back("GENIE Background");
      SampleTypes.push_back("Background");
      SampleFiles.push_back("run1_FHC/analysisOutputFHC_Overlay_GENIE_Background_All.root");

      SampleNames.push_back("GENIE Hyperon");
      SampleTypes.push_back("Hyperon");
      SampleFiles.push_back("run1_FHC/analysisOutputFHC_Overlay_GENIE_Hyperon_All.root");

      SelectionParameters P = P_FHC_Tune_325_NoBDT;

      EventAssembler E;
      SelectionManager M(P);
      M.SetPOT(POT);
      M.ImportSelectorBDTWeights(P.p_SelectorBDT_WeightsDir);

      // Setup the histograms
      M.SetupHistograms(100,0.0,5.0,";Neutrino Energy (GeV);Events");

      // Sample Loop
      for(size_t i_s=0;i_s<SampleNames.size();i_s++){

         E.SetFile(SampleFiles.at(i_s));
         if(SampleTypes.at(i_s) != "EXT" && SampleTypes.at(i_s) != "Data") M.AddSample(SampleNames.at(i_s),SampleTypes.at(i_s),E.GetPOT());
         else if(SampleTypes.at(i_s) == "Data") M.AddSample(SampleNames.at(i_s),SampleTypes.at(i_s),Data_POT);
         else if(SampleTypes.at(i_s) == "EXT") M.AddSample(SampleNames.at(i_s),SampleTypes.at(i_s),EXT_POT);

         // Event Loop
         for(int i=0;i<E.GetNEvents();i++){

            if(i % 10000 == 0) std::cout << i << "/" << E.GetNEvents() << std::endl;

            Event e = E.GetEvent(i);

            M.SetSignal(e);                
            M.AddEvent(e);

            if(!M.FiducialVolumeCut(e)) continue;
            if(!M.TrackCut(e)) continue;
            if(!M.ShowerCut(e)) continue;
            if(!M.ChooseMuonCandidate(e)) continue;
            //if(!M.ChooseProtonPionCandidates(e)) continue;

            double E = e.Neutrino.at(0).E;

            M.FillHistograms(e,E);                
         }
         E.Close();
      }

      M.DrawHistograms(label);
   }
