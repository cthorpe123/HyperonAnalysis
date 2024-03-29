#ifndef _SidebandFitter_cxx_
#define _SidebandFitter_cxx_

#include "SidebandFitter.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// Constructor

SidebandFitter::SidebandFitter()
{

   Minimizer->SetMaxFunctionCalls(100);
   Minimizer->SetTolerance(0.01);

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Add the histograms with MC and data to be used in the fit
// h_fixed_v = vector of histograms to be kept fixed in the fit
// h_tovar_v = vector of histograms to be varied in the fit to find optimal scaling on each
// h_data = data to be used in the fit

void SidebandFitter::AddHistograms(std::vector<TH1D*> h_fixed_v,std::vector<TH1D*> h_tovar_v,TH1D* h_data){

   AddFixedHistograms(h_fixed_v);
   AddVariedHistograms(h_tovar_v);
   AddDataHistogram(h_data);

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Add the histograms to be kept fixed in the fit

void SidebandFitter::AddFixedHistograms(std::vector<TH1D*> h_fixed_v){

   for(int i=0;i<h_fixed_v.size();i++){   
      h_Fixed_v.push_back((TH1D*)h_fixed_v.at(i)->Clone(Form("h_Fixed_%",i)));
      h_Fixed_v.back()->SetDirectory(0);
   }

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Add the histograms to be varied in the fit to find optimal scaling on each

void SidebandFitter::AddVariedHistograms(std::vector<TH1D*> h_tovar_v){

   for(int i=0;i<h_tovar_v.size();i++){   
      h_ToVar_v.push_back((TH1D*)h_tovar_v.at(i)->Clone(Form("h_ToVar_%",i)));
      h_ToVar_v.back()->SetDirectory(0);
   }

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Add the histogram containing the data to be used in the fit

void SidebandFitter::AddDataHistogram(TH1D* h_data){

   h_Data = (TH1D*)h_data->Clone("h_Data");
   h_Data->SetDirectory(0);

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Perform the fit

FitResult SidebandFitter::DoFit(){

   if(h_Data == nullptr) 
      throw std::invalid_argument("SidebandFitter: Data histogram not loaded!");
   if(!h_ToVar_v.size()) 
      throw std::invalid_argument("SidebandFitter: No histograms set to be varied in fit");

   const int nbins = h_Data->GetNbinsX();
   FitResult result;

   // Calculate the total contribution from the non-varied histograms
   TH1D* h_tot_fixed = (TH1D*)h_Data->Clone("h_tot_fixed");
   h_tot_fixed->Reset();
   for(int i=0;i<h_Fixed_v.size();i++) h_tot_fixed->Add(h_Fixed_v.at(i));     
/*
   TH1D* h_tot_varied = (TH1D*)h_ToVar_v.at(0)->Clone("h_tot_varied");
   for(size_t i=1;i<h_ToVar_v.size();i++) h_tot_varied->Add(h_ToVar_v.at(i));
   std::cout << "Before Fit" << std::endl;
   for(int i=1;i<h_Data->GetNbinsX()+1;i++)
      std::cout << "Data = " << h_Data->GetBinContent(i) << "  Fixed = " << h_tot_fixed->GetBinContent(i) << " Varied = " << h_tot_varied->GetBinContent(i) << std::endl; 
  */

   // Set function
   ROOT::Math::Functor min = ROOT::Math::Functor([&](const double *coeff){

         TMatrixDSym covstat(nbins);
         for(int i_b=1;i_b<nbins+1;i_b++){
         covstat[i_b-1][i_b-1] += h_Data->GetBinError(i_b)*h_Data->GetBinError(i_b) + h_tot_fixed->GetBinError(i_b)*h_tot_fixed->GetBinError(i_b);
         for(int i_h=0;i_h<h_ToVar_v.size();i_h++)
         covstat[i_b-1][i_b-1] += h_ToVar_v.at(i_h)->GetBinError(i_b)*h_ToVar_v.at(i_h)->GetBinError(i_b)*coeff[i_h]*coeff[i_h];
         }

         TMatrixDSym covinv = covstat;
         covinv.Invert();

         double chi2 = 0;
         double ndof=nbins-1;
         for(int i=1;i<nbins+1;i++){         
         for(int j=1;j<nbins+1;j++){         
         double diff_i = h_Data->GetBinContent(i) - h_tot_fixed->GetBinContent(i);
         double diff_j = h_Data->GetBinContent(j) - h_tot_fixed->GetBinContent(j);        
         for(int i_h=0;i_h<h_ToVar_v.size();i_h++){
         diff_i -= h_ToVar_v.at(i_h)->GetBinContent(i)*coeff[i_h];
         diff_j -= h_ToVar_v.at(i_h)->GetBinContent(j)*coeff[i_h];
         }   
         chi2 += diff_i*covinv[i-1][j-1]*diff_j;
         }
         }

         return chi2/ndof;
   },h_ToVar_v.size());    

   Minimizer->SetFunction(min);
   for(size_t i=0;i<h_ToVar_v.size();i++){
      std::string name = "scale_" + std::to_string(i);
      Minimizer->SetVariable(i,name,1.0,2.0);
   }

   int FitStatus = Minimizer->Minimize(); 
   result.FitCov.ResizeTo((Int_t)h_ToVar_v.size(),(Int_t)h_ToVar_v.size());

   for(size_t i=0;i<h_ToVar_v.size();i++){
      result.FitVals.push_back(Minimizer->X()[i]);
      for(size_t j=0;j<h_ToVar_v.size();j++)
         result.FitCov[i][j] = Minimizer->CovMatrix(i,j);
   }

   return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Reset the fitter (not the same as destructor)

void SidebandFitter::Clear(){

   h_ToVar_v.clear();
   h_Fixed_v.clear();
   delete h_Data;

}

///////////////////////////////////////////////////////////////////////////////////////////////

#endif
