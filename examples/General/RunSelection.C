R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so)
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so)

#include "SelectionManager.h"
#include "EventAssembler.h"
#include "Cut.h"
#include "SelectionParameters.h"

#include "Parameters.h"
#include "SampleSets_Example.h"

   void RunSelection(){

      double POT = 1.0e21; // POT to scale samples to

      BuildTunes();
      SelectionParameters P = P_FHC_Tune_325_NoBDT;

      std::string label = "test";

      // Setup selection manager. Set POT to scale sample to, import the BDT weights
      EventAssembler E;
      SelectionManager M(P);
      M.SetPOT(POT);
      M.ImportSelectorBDTWeights(P.p_SelectorBDT_WeightsDir);
      //M.ImportAnalysisBDTWeights(P.p_AnalysisBDT_WeightsDir);

      double Signal = 0.0;
      double BG = 0.0;
      double Sel_Signal = 0.0;
      double Sel_BG = 0.0;

      E.SetFile("run1_FHC/analysisOutputFHC_Overlay_GENIE_Hyperon_All.root");
      M.AddSample("Hyperon","Hyperon",E.GetPOT());

      // Event Loop
      for(int i=0;i<E.GetNEvents();i++){

         if(i % 10000 == 0) std::cout << i << "/" << E.GetNEvents() << std::endl;

         Event e = E.GetEvent(i);

         M.SetSignal(e);                
         M.AddEvent(e);

         if(e.EventIsSignal) Signal += e.Weight;
         else BG += e.Weight;

         if(!M.FiducialVolumeCut(e)) continue;
         if(!M.TrackCut(e)) continue;
         if(!M.ShowerCut(e)) continue;
         if(!M.ChooseMuonCandidate(e)) continue;
         if(!M.ChooseProtonPionCandidates(e)) continue;
         //if(!M.AnalysisBDTCut(e)) continue;       
         if(!M.ConnectednessTest(e)) continue;           
         if(!M.WCut(e)) continue;
         if(!M.AngleCut(e)) continue;

         if(e.EventIsSignal) Sel_Signal += e.Weight;
         else Sel_BG += e.Weight;
      }

      E.Close();

      std::cout << "Signal = " << Signal << "  Selected = " << Sel_Signal << std::endl;
      std::cout << "BG = " << BG << "  Selected = " << Sel_BG << std::endl;
   }
