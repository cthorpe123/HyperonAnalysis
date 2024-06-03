#ifndef _PlottingFunctions_h_
#define _PlottingFunctions_h_

#include "TFile.h"
#include "TH1D.h"
#include "Misc.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "FluxWeight2.h"
#include "TMatrixDSym.h"

using std::string;
using std::vector;

const double _EPSILON_ = 1e-5;

namespace HypPlot {

   bool DrawWatermark = true;

   // Single panel axis settings etc.

   const double Single_CanvasX = 800;
   const double Single_CanvasY = 600; 
   const double Single_PadSplit = 0.85; // Height at which to split the cavas into plot/legend

   const double Single_XaxisTitleSize = 0.05;
   const double Single_YaxisTitleSize = 0.05;
   const double Single_XaxisTitleOffset = 0.93;
   const double Single_YaxisTitleOffset = 1.06;
   const double Single_XaxisLabelSize = 0.045;
   const double Single_YaxisLabelSize = 0.045;

   const double Single_TextLabelSize = 0.09; // Label size if bin labels are used

   // Two panel axis settings etc.

   const double Dual_CanvasX = 800;
   const double Dual_CanvasY = 750; 
   const double Dual_PadSplitLow = 0.3; // Height at which to split the canvas between the ratio and the main plot
   const double Dual_PadSplitHigh = 0.9; // Height at which to split the canvas between the main plot and the legend

   const double Dual_MainXaxisTitleSize = 0;
   const double Dual_MainYaxisTitleSize = 0.065;
   const double Dual_MainXaxisTitleOffset = 0;
   const double Dual_MainYaxisTitleOffset = 0.8;
   const double Dual_MainXaxisLabelSize = 0;
   const double Dual_MainYaxisLabelSize = 0.05;

   const double Dual_RatioXaxisTitleSize = 0.12;
   const double Dual_RatioYaxisTitleSize = 0.12;
   const double Dual_RatioXaxisTitleOffset = 0.8;
   const double Dual_RatioYaxisTitleOffset = 0.43;
   const double Dual_RatioXaxisLabelSize = 0.1;
   const double Dual_RatioYaxisLabelSize = 0.1;

   const double Dual_TextLabelSize = 0.17; // Label size if bin labels are used

   // Matrix axis settings etc.

   const double Matrix_CanvasX = 800;
   const double Matrix_CanvasY = 600; 

   const double Matrix_XaxisTitleSize = 0.05;
   const double Matrix_YaxisTitleSize = 0.05;
   const double Matrix_XaxisTitleOffset = 0.93;
   const double Matrix_YaxisTitleOffset = 1.02;
   const double Matrix_XaxisLabelSize = 0.045;
   const double Matrix_YaxisLabelSize = 0.045;
   const double Matrix_ZaxisLabelSize = 0.045;

   const double Matrix_TextLabelSize = 0.07;

   const int GoodLineColors[13] = {1,2,8,4,6,38,46,43,30,9,7,14,3};

///////////////////////////////////////////////////////////////////////////////////////////////

void DontDrawWatermark(){

   DrawWatermark = false;

}

///////////////////////////////////////////////////////////////////////////////////////////////

// Calculate the chi2 score comparing h_Pred to h_Data with cov encoding any systematic
// uncertainties.
// Inputs:
// h_Pred = 1D hist holding predictions
// h_Data = 1D hist with data
// cov = covariance matrix with systematic uncertainties only
// skip = any bins not to be included in calculation 

std::pair<double,int> Chi2(TH1D* h_Pred,TH1D* h_Data,TMatrixDSym cov=TMatrixDSym(0),std::vector<int> skip={}){
 
   const int nbins = h_Data->GetNbinsX();

   // Remove any bins with no events in
   std::vector<int> nonzero_bins;
   for(size_t i=1;i<nbins+1;i++)
      if(h_Pred->GetBinContent(i) > 0 && std::find(skip.begin(),skip.end(),i) == skip.end()) nonzero_bins.push_back(i);

   TMatrixDSym cov_nonzero(nonzero_bins.size());

   if(cov.GetNcols()){
      if(cov.GetNcols() != nbins) 
         throw std::invalid_argument("PlottingFunctions: Attempting to get get chi2 score when covariance matrix has size " + std::to_string(cov.GetNcols()) + "x" + std::to_string(cov.GetNcols()) + " hists have size "  + std::to_string(nbins));
      for(size_t i=0;i<nonzero_bins.size();i++)
         for(size_t j=0;j<nonzero_bins.size();j++)
            cov_nonzero[i][j] = cov[nonzero_bins.at(i)-1][nonzero_bins.at(j)-1];        
   }

   // Add the statistical uncertainty
   TMatrixDSym cov_stat(nonzero_bins.size());
   for(size_t i=0;i<nonzero_bins.size();i++)
      cov_stat[i][i] = h_Pred->GetBinError(nonzero_bins.at(i))*h_Pred->GetBinError(nonzero_bins.at(i)) + h_Data->GetBinError(nonzero_bins.at(i))*h_Data->GetBinError(nonzero_bins.at(i));

   cov_nonzero += cov_stat;
   cov_nonzero.Invert();

   double Chi2 = 0.0;
   for(size_t i=0;i<nonzero_bins.size();i++){
      for(size_t j=0;j<nonzero_bins.size();j++){
         Chi2 += (h_Pred->GetBinContent(nonzero_bins.at(i)) - h_Data->GetBinContent(nonzero_bins.at(i)))*cov_nonzero[i][j]*(h_Pred->GetBinContent(nonzero_bins.at(j)) - h_Data->GetBinContent(nonzero_bins.at(j)));              
}
}

   std::cout << "Chi2/ndof = " << Chi2 << "/" << nonzero_bins.size() << std::endl;

   return {Chi2,nonzero_bins.size()};
}

///////////////////////////////////////////////////////////////////////////////////////////////
/*
// Calculate the chi2 score with a prediction, data, and covariance matrix describing any
// systematic uncertainties

std::pair<double,int> Chi2(const TH1D* h_pred,const TH1D* h_data,TMatrixDSym* cov=nullptr){

   const int nbins = h_pred->GetNbinsX();
   if(cov == nullptr) cov = new TMatrixDSym(nbins);

   TMatrixDSym cov_stat(nbins);
   for(int i=1;i<nbins+1;i++)
      cov_stat[i-1][i-1] = h_pred->GetBinError(i)*h_pred->GetBinError(i) + h_data->GetBinError(i)*h_data->GetBinError(i);
   cov_stat += *cov;
   cov_stat.Invert();

   double chi2 = 0.0;
   for(int i=1;i<nbins+1;i++){
      for(int j=1;j<nbins+1;j++){
         double diff_i = h_pred->GetBinContent(i) - h_data->GetBinContent(i);
         double diff_j = h_pred->GetBinContent(j) - h_data->GetBinContent(j);
         chi2 += diff_i*diff_j*cov_stat[i-1][j-1];
      }
   }

   return std::make_pair(chi2,nbins-1);
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////

// Take set of histograms to be stacked together and compute the total stat error
// h_Cov = 2D histogram storing the covariance matrix, systematic errors only

TH1D * MakeErrorBand(std::vector<TH1D*> hists,TH2D* h_Cov=nullptr){

   TH1D *h_errors = (TH1D*)hists.at(0)->Clone("h_errors");  
   h_errors->Reset();

   // Iterate over the bins of the contributing histograms
   for(int i_b=1;i_b<h_errors->GetNbinsX()+1;i_b++){
      double events = 0.0;
      double variance = 0.0;
      for(size_t i_h=0;i_h<hists.size();i_h++){
         events += hists.at(i_h)->GetBinContent(i_b);
         variance += hists.at(i_h)->GetBinError(i_b)*hists.at(i_h)->GetBinError(i_b);
      }
      if(h_Cov != nullptr) variance += h_Cov->GetBinContent(i_b,i_b);
      h_errors->SetBinContent(i_b,events);
      h_errors->SetBinError(i_b,sqrt(variance));
   }

   return h_errors;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the maximum of a histogram including error bars

double GetHistMaxError(TH1D *h){

   double max = -1e10;

   for(int i_b=1;i_b<h->GetNbinsX()+1;i_b++)
      if(h->GetBinContent(i_b)+h->GetBinError(i_b) > max)
         max = h->GetBinContent(i_b)+h->GetBinError(i_b);

   return max;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Get the maximum of a histogram not including the errors

double GetHistMax(TH1D *h){

   double max = -1e10;

   for(int i_b=1;i_b<h->GetNbinsX()+1;i_b++)
      if(h->GetBinContent(i_b) > max)
         max = h->GetBinContent(i_b);

   return max;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Make a 2D histogram from a TMatrix object
// Uses h_example for the bin widths, axis titles etc.

TH2D* MakeHistogram(TMatrixD Mat,TH2D* h_example,string Name=""){

   TH2D* h_Mat = (TH2D*)h_example->Clone(("h_" + Name).c_str());
   const int nbins = h_Mat->GetNbinsX();

   for(int i=1;i<nbins+1;i++)
      for(int j=1;j<nbins+1;j++)
         h_Mat->SetBinContent(i,j,Mat[i-1][j-1]);

   return h_Mat; 
}

///////////////////////////////////////////////////////////////////////////////////////////////
/*
void DrawMatrix(TMatrixD Mat,TH2D* h_example,string title,bool uselabels=false,bool usetext=false){

   TH2D* h = MakeHistogram(Mat,h_example);
   DrawMatrix(h,h_example,title,uselabels,usetext);

}
*/
///////////////////////////////////////////////////////////////////////////////////////////////
// Make a TGraph object from a vector containing the x positions of data events in 
// some variable

TGraph* MakeDataGraph(std::vector<double> data_v){

   std::vector<Double_t> Y(data_v.size(),0.0); 

   TGraph* g = new TGraph(data_v.size(),&(data_v[0]),&(Y[0]));
   g->SetMarkerStyle(23);
   g->SetMarkerSize(3);

   return g;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Draw a 2D matrix as a histogram
// Uses h_example for the bin widths, axis titles etc.
// Can import bin labels from h_example and write bin content on bins if corresponding bools
// are set.

void DrawMatrix(TH2D* h,TH2D* h_example,string title,std::string plotdir,bool uselabels=false,bool usetext=false){

   h->SetContour(1000);
   h->GetXaxis()->SetTitleSize(Matrix_XaxisTitleSize);
   h->GetYaxis()->SetTitleSize(Matrix_YaxisTitleSize);
   h->GetXaxis()->SetTitleOffset(Matrix_XaxisTitleOffset);
   h->GetYaxis()->SetTitleOffset(Matrix_YaxisTitleOffset);
   h->GetXaxis()->SetLabelSize(Matrix_XaxisLabelSize);
   h->GetYaxis()->SetLabelSize(Matrix_YaxisLabelSize);
   h->GetZaxis()->SetLabelSize(Matrix_ZaxisLabelSize);

   TCanvas *c = new TCanvas("c","c",Matrix_CanvasX,Matrix_CanvasY);
   c->SetRightMargin(0.13);

   if(uselabels){
      c->SetLeftMargin(0.15);
      h->GetXaxis()->SetLabelSize(Matrix_TextLabelSize);
      h->GetYaxis()->SetLabelSize(Matrix_TextLabelSize);
   }

   // Generates text with bin content, uses 3 sf
   if(usetext){
      const int sf = 3;
      c->SetRightMargin(0.05);
      h->Draw("col");
      h->SetMarkerSize(3);
      for(int i=1;i<h->GetNbinsX()+1;i++){
         for(int j=1;j<h->GetNbinsX()+1;j++){
            double content = setsf(h->GetBinContent(i,j),sf); 
            string text = std::to_string(content);        
            if(abs(content) < 1e-4 || abs(content) > 1e4){
               int OM = floor(log10(abs(content)));               
               text = to_string_with_precision(content*pow(10,-OM),2) + " #times 10^{" + std::to_string(OM) + "}";
            }
            if(abs(content) < 1e-10) text = std::to_string(0);
            auto t = new TLatex(h->GetXaxis()->GetBinCenter(i),h->GetYaxis()->GetBinCenter(j),text.c_str());
            t->SetTextAlign(22);
            t->Draw();
            t->SetTextSize(0.04);
         }
      }
   }
   else h->Draw("colz");

   TLegend *l_Watermark;
   if(!usetext) l_Watermark = new TLegend(0.445,0.91,0.910,0.99);
   else l_Watermark = new TLegend(0.485,0.91,0.950,0.99);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.001);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextFont(62);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");

   h->SetStats(0);
   if(DrawWatermark) l_Watermark->Draw();

   system(("mkdir -p " + plotdir).c_str());   
   c->Print((plotdir + "/" + title + ".png").c_str());
   c->Print((plotdir + "/" + title + ".pdf").c_str());
   c->Print((plotdir + "/" + title + ".C").c_str());
   c->Close();

   //delete h;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Draw a stack of histograms
// If you anto to include systematic uncertainties in the MC, make sure they're already included
// in the error of h_errors
// Inputs:
// hist_v = the stack of histograms you want to draw
// h_errors = the total MC prediction with uncertainties
// captions = list of what should be added to the legend for each histogram
// plotdir = the directory to write the plots into
// label = label to be added the front of the plot
// mode = FHC/RHC/BNB
// run = 1,3 or 13 (for 1+3)
// POT = the POT used 
// signalscale = scaling for the signal
// hasdata = flag indicating if there is data in the stack of histograms
// colors = list of colors for the histograms
// binlabels = bin labels
// chi2ndof = chi2 and degrees of freedom if applicable
// data_v list of x positions of individual data events (used for signal box plots only)

void DrawHistogram(std::vector<TH1D*> hist_v,TH1D* h_errors,TH1D* h_data,vector<string> captions,string plotdir,string label,vector<int> mode,vector<int> run,vector<double> POT,double signalscale,bool hasdata,vector<int> colors,vector<string> binlabels,std::pair<double,int> chi2ndof,std::vector<double> data_v={}){

   bool drawchi2 = chi2ndof.second != 0;
 
   assert(mode.size() == run.size() && run.size() == POT.size() && mode.size() < 3);   
   for(size_t i_r=0;i_r<run.size();i_r++) assert(mode.at(i_r) == kFHC || mode.at(i_r) == kRHC || mode.at(i_r) == kBNB);
   if(hasdata && h_data == nullptr) throw std::invalid_argument("PlottingFunctions::DrawHistogram: hasdata flag set to true but data histogram is nullptr, exiting");

   const int nbins = hist_v.at(0)->GetNbinsX();

   // Set the bin labels
   if(binlabels.size()){
      for(int i=1;i<nbins+1;i++){
         for(size_t i_h=0;i_h<hist_v.size();i_h++) 
            hist_v.at(i_h)->GetXaxis()->SetBinLabel(i,binlabels.at(i-1).c_str());
         h_errors->GetXaxis()->SetBinLabel(i,binlabels.at(i-1).c_str());
         if(hasdata) h_data->GetXaxis()->SetBinLabel(i,binlabels.at(i-1).c_str());
      }
   }

   // Setup the canvas
   TCanvas *c = new TCanvas("c","c",Single_CanvasX,Single_CanvasY);
   TPad *p_plot = new TPad("p_plot","p_plot",0,0,1,Single_PadSplit);
   TPad *p_legend = new TPad("p_legend","p_legend",0,Single_PadSplit,1,1);
   p_legend->SetBottomMargin(0);
   p_legend->SetTopMargin(0.1);
   p_plot->SetTopMargin(0.01);

   THStack *hs = new THStack("hs",hist_v.at(0)->GetTitle());

   // Create the empty legend
   TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
   l->SetBorderSize(0);
   const int nhists = hist_v.size() + static_cast<int>(hasdata);
   int ncols = 2;
   if(nhists > 4) ncols = 3;
   if(nhists > 6) ncols = 4;
   if(nhists > 12) ncols = 5;
   l->SetNColumns(ncols);

   for(size_t i_h=0;i_h<hist_v.size();i_h++){ 
      hist_v.at(i_h)->SetFillColor(colors.at(i_h));
      hs->Add(hist_v.at(i_h),"HIST");
      if(hasdata) l->AddEntry(hist_v.at(i_h),(captions.at(i_h) + " = " + to_string_with_precision(hist_v.at(i_h)->Integral(),1)).c_str(),"F");
      else l->AddEntry(hist_v.at(i_h),captions.at(i_h).c_str(),"F");
   }

   if(hasdata){
      h_data->SetLineWidth(1);
      h_data->SetLineColor(1);
      h_data->SetMarkerStyle(20);
      h_data->SetMarkerColor(1);
      l->AddEntry(h_data,("Data = "+ to_string_with_precision(h_data->Integral(),1)).c_str(),"P");
   }

   // Create the "MicroBooNE" watermark
   TLegend *l_Watermark = new TLegend(0.45,0.900,0.89,0.985);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.005);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetTextFont(62);

   if(!hasdata && !data_v.size()) l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");
   else if(run.size() == 1) l_Watermark->SetHeader(("MicroBooNE Run " + std::to_string(run.at(0)) + ", Preliminary").c_str());
   else if(run.size() == 2) l_Watermark->SetHeader(("MicroBooNE Runs " + std::to_string(run.at(0)) + " + " + std::to_string(run.at(1)) + ", Preliminary").c_str());
   else throw std::invalid_argument("PlottingFunctions::DrawHistogram: Currently maximum of two running periods supported");
/*
   if(!hasdata && !data_v.size()) l_Watermark->SetHeader("MicroBooNE Simulation","R");
   else if(run.size() == 1) l_Watermark->SetHeader(("MicroBooNE Run " + std::to_string(run.at(0))).c_str());
   else if(run.size() == 2) l_Watermark->SetHeader(("MicroBooNE Runs " + std::to_string(run.at(0)) + " + " + std::to_string(run.at(1))).c_str());
   else throw std::invalid_argument("PlottingFunctions::DrawHistogram: Currently maximum of two running periods supported");
*/

   // Create the POT label
   //TLegend *l_POT = new TLegend(0.54,0.815,0.885,0.900);
   TLegend *l_POT,*l_POT2;
   if(DrawWatermark) l_POT = new TLegend(0.54,0.815,0.885,0.900);
   else l_POT = new TLegend(0.54,0.900,0.885,0.985);
   l_POT->SetBorderSize(0);
   l_POT->SetMargin(0.005);
   l_POT->SetTextAlign(32);
   if(DrawWatermark) l_POT2 = new TLegend(0.536,0.735,0.886,0.815,NULL,"brNDC");
   else l_POT2 = new TLegend(0.536,0.815,0.886,0.900,NULL,"brNDC");
   //TLegend *l_POT2 = new TLegend(0.54,0.735,0.89,0.815,NULL,"brNDC");
   l_POT2->SetBorderSize(0);
   l_POT2->SetMargin(0.005);
   l_POT2->SetTextAlign(32);

   // Chi2
   TLegend *l_Chi2 = new TLegend(0.12,0.900,0.36,0.985);
   l_Chi2->SetBorderSize(0);
   l_Chi2->SetMargin(0.005);
   l_Chi2->SetTextAlign(12);
   l_Chi2->SetTextSize(0.05);
   l_Chi2->SetHeader(("#chi^{2}/ndof = " + to_string_with_precision(chi2ndof.first,1) + "/" + std::to_string(chi2ndof.second)).c_str());

   if(mode.at(0) == kFHC) l_POT->SetHeader(("NuMI FHC, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.at(0) == kRHC) l_POT->SetHeader(("NuMI RHC, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.at(0) == kBNB) l_POT->SetHeader(("BNB, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());

   if(mode.size() == 2 && mode.at(1) == kFHC) l_POT2->SetHeader(("NuMI FHC, " + to_string_with_precision(POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.size() == 2 && mode.at(1) == kRHC) l_POT2->SetHeader(("NuMI RHC, " + to_string_with_precision(POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.size() == 2 && mode.at(0) == kBNB && mode.at(1) == kBNB) l_POT->SetHeader(("BNB, " + to_string_with_precision(POT.at(0) + POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());

   // Setup the axis titles etc
   h_errors->SetTitle(hs->GetTitle());
   h_errors->SetFillStyle(3253);
   h_errors->SetFillColor(1);
   h_errors->GetXaxis()->SetTitleSize(Single_XaxisTitleSize);
   h_errors->GetYaxis()->SetTitleSize(Single_YaxisTitleSize);
   h_errors->GetXaxis()->SetTitleOffset(Single_XaxisTitleOffset);
   h_errors->GetYaxis()->SetTitleOffset(Single_YaxisTitleOffset);
   h_errors->GetXaxis()->SetLabelSize(Single_XaxisLabelSize);
   h_errors->GetYaxis()->SetLabelSize(Single_YaxisLabelSize);
   if(binlabels.size()) h_errors->GetXaxis()->SetLabelSize(Single_TextLabelSize);

   // Draw everything
   p_legend->Draw();
   p_legend->cd();
   l->Draw();
   c->cd();
   p_plot->Draw();
   p_plot->cd();

   h_errors->Draw("E2");
   hs->Draw("HIST same");
   h_errors->Draw("E2 same");
   if(hasdata) h_data->Draw("E0 P0 same");
   double maximum = GetHistMaxError(h_errors);
   if(hasdata) maximum = std::max(maximum,GetHistMaxError(h_data));
   h_errors->GetYaxis()->SetRangeUser(0.0,maximum*1.25);
   h_errors->SetStats(0);

   // Draw the data graph if required
   TGraph* g_data1;
   TGraph* g_data2;
   if(data_v.size()){
      std::vector<double> c_height1(data_v.size(),maximum*0.039);
      std::vector<double> c_height2(data_v.size(),maximum*0.075);
      TGraph* g_data1 = new TGraph(data_v.size(),&(data_v[0]),&(c_height1[0]));
      TGraph* g_data2 = new TGraph(data_v.size(),&(data_v[0]),&(c_height2[0]));
      g_data1->SetMarkerStyle(23);
      g_data1->SetMarkerSize(3);
      g_data2->SetMarkerStyle(2);
      g_data2->SetMarkerSize(3);
      g_data1->Draw("P same");   
      g_data2->Draw("P same");   
      l->AddEntry(g_data1,"Events","P"); 
   }

   p_plot->RedrawAxis();

   // Draw the various legends, labels etc.
   if(POT.size() > 0) l_POT->Draw();
   if(POT.size() == 2 && mode.at(0) == kFHC && mode.at(1) == kRHC){
      l_POT2->Draw();
      h_errors->GetYaxis()->SetRangeUser(0.0,maximum*1.4);
   }
   if(DrawWatermark) l_Watermark->Draw();
   if(hasdata && drawchi2) l_Chi2->Draw();

   c->cd();
   system(("mkdir -p " + plotdir).c_str());
   c->Print((plotdir + "/" + label + ".png").c_str());
   c->Print((plotdir + "/" + label + ".pdf").c_str());
   c->Print((plotdir + "/" + label + ".C").c_str());
   c->Clear();
   c->Close();

   // Draw the ratio and get chi2 if there is data
   if(hasdata){

      // Setup split canvas
      TCanvas *c2 = new TCanvas("c2","c2",Dual_CanvasX,Dual_CanvasY);
      TPad *p_ratio2 = new TPad("p_ratio2","p_ratio2",0,0.0,1,Dual_PadSplitLow);
      TPad *p_plot2 = new TPad("p_plot2","p_plot2",0,Dual_PadSplitLow,1,Dual_PadSplitHigh);
      TPad *p_legend2 = new TPad("p_legend2","p_plot2",0,Dual_PadSplitHigh,1,1);
      p_legend2->SetBottomMargin(0);
      p_legend2->SetTopMargin(0.1);
      p_ratio2->SetTopMargin(0.005);
      p_plot2->SetTopMargin(0.01);
      p_plot2->SetBottomMargin(0.02);
      p_ratio2->SetGrid(0,1);
      p_ratio2->SetBottomMargin(0.22);
      p_legend2->Draw();
      p_legend2->cd();
      l->Draw();
      c2->cd();
      p_plot2->Draw();
      p_plot2->cd();

      h_errors->GetXaxis()->SetLabelSize(Dual_MainXaxisLabelSize);
      h_errors->GetYaxis()->SetLabelSize(Dual_MainYaxisLabelSize);
      h_errors->GetXaxis()->SetTitleSize(Dual_MainXaxisTitleSize);
      h_errors->GetYaxis()->SetTitleSize(Dual_MainYaxisTitleSize);
      h_errors->GetYaxis()->SetTitleOffset(Dual_MainYaxisTitleOffset);

      // Draw the main plot
      h_errors->Draw("E2");
      hs->Draw("HIST same");    
      h_errors->Draw("E2 same");
      h_data->Draw("e0 same");

      // Draw the various legends, labels etc.
      if(POT.size() > 0) l_POT->Draw();
      if(POT.size() == 2 && mode.at(0) == kFHC && mode.at(1) == kRHC) l_POT2->Draw();
      if(DrawWatermark) l_Watermark->Draw();
      if(drawchi2) l_Chi2->Draw();

      // Make the raio plot
      TH1D *h_errors_ratio = (TH1D*)h_errors->Clone("h_errors_ratio");
      TH1D* h_data_ratio = (TH1D*)h_data->Clone("h_data_ratio");

      std::pair<double,double> minmax = {1.0,1.0};

      for(int i=1;i<h_errors_ratio->GetNbinsX()+1;i++){

         double pred = h_errors->GetBinContent(i);
         double mc_error = h_errors->GetBinError(i)/pred;
         double data_ratio = h_data->GetBinContent(i)/pred;
         double data_ratio_error = h_data->GetBinError(i)/pred;

         h_errors_ratio->SetBinContent(i,1.0);

         if(pred > 0){
            h_errors_ratio->SetBinError(i,mc_error);            
            h_data_ratio->SetBinContent(i,data_ratio);
            h_data_ratio->SetBinError(i,data_ratio_error);            
            minmax.first = std::min(minmax.first,std::min(1-mc_error,data_ratio-data_ratio_error));
            minmax.second = std::max(minmax.second,std::max(1+mc_error,data_ratio+data_ratio_error));                          
         }
         else{
            h_errors_ratio->SetBinError(i,0.0); 
            h_data_ratio->SetBinContent(i,1.0); 
            h_data_ratio->SetBinError(i,0.0); 
         }
      }

      // Draw the ratio plot
      c2->cd();
      p_ratio2->Draw();
      p_ratio2->cd();
      h_errors_ratio->Draw("E2");
      h_errors_ratio->GetYaxis()->SetRangeUser(0.95*minmax.first,1.05*minmax.second);

      h_errors_ratio->SetStats(0);
      //h_errors_ratio->GetXaxis()->SetTitle(hist_v.at(0)->GetXaxis()->GetTitle());
      h_errors_ratio->GetYaxis()->SetTitle("Data/MC");

      h_errors_ratio->GetXaxis()->SetTitleSize(Dual_RatioXaxisTitleSize);
      h_errors_ratio->GetYaxis()->SetTitleSize(Dual_RatioYaxisTitleSize);
      h_errors_ratio->GetXaxis()->SetTitleOffset(Dual_RatioXaxisTitleOffset);
      h_errors_ratio->GetYaxis()->SetTitleOffset(Dual_RatioYaxisTitleOffset);
      h_errors_ratio->GetXaxis()->SetLabelSize(Dual_RatioXaxisLabelSize);
      h_errors_ratio->GetYaxis()->SetLabelSize(Dual_RatioYaxisLabelSize);
      if(binlabels.size()){
         h_errors_ratio->GetXaxis()->SetLabelOffset(0.02);
         h_errors_ratio->GetXaxis()->SetLabelSize(Dual_TextLabelSize);
      }

      h_data_ratio->Draw("E0 same");
      p_plot2->RedrawAxis();

      c2->Print((plotdir + "/" + label + "_Ratio.pdf").c_str());
      c2->Print((plotdir + "/" + label + "_Ratio.png").c_str());
      c2->Print((plotdir + "/" + label + "_Ratio.C").c_str());
      c2->Close();
   }

   c->Close();

}

///////////////////////////////////////////////////////////////////////////////////////////////

// Draw a stack of histograms
// If you anto to include systematic uncertainties in the MC, make sure they're already included
// in the error of h_errors
// Inputs:
// hist_v = the stack of histograms you want to draw (including data if applicable)
// h_errors = the total MC prediction with uncertainties
// captions = list of what should be added to the legend for each histogram
// plotdir = the directory to write the plots into
// label = label to be added the front of the plot
// mode = FHC/RHC/BNB
// run = 1,3 or 13 (for 1+3)
// POT = the POT used 
// signalscale = scaling for the signal
// binlabels = bin labels

void DrawHistogramSys(vector<TH1D*> hist_v,TH1D* h_cv,string plotdir,string label,string type,string dialname,vector<int> mode,vector<int> run,vector<double> POT,vector<string> binlabels){

   assert(mode.size() == run.size() && run.size() == POT.size() && mode.size() < 3);   
   for(size_t i_r=0;i_r<run.size();i_r++) assert(mode.at(i_r) == kFHC || mode.at(i_r) == kRHC || mode.at(i_r) == kBNB);

   // Set the bin labels
   if(binlabels.size())
      for(size_t i_u=0;i_u<hist_v.size();i_u++)
         for(int i=1;i<h_cv->GetNbinsX()+1;i++) hist_v.at(i_u)->GetXaxis()->SetBinLabel(i,binlabels.at(i-1).c_str());

   // Setup the canvas and pads
   TCanvas *c = new TCanvas("c","c",Single_CanvasX,Single_CanvasY);
   TPad *p_plot = new TPad("pad1","pad1",0,0,1,Single_PadSplit);
   TPad *p_legend = new TPad("pad2","pad2",0,Single_PadSplit,1,1);
   p_legend->SetBottomMargin(0);
   p_legend->SetTopMargin(0.1);
   p_plot->SetTopMargin(0.01);

   std::string title = ";" + string(h_cv->GetXaxis()->GetTitle()) + ";" + string(h_cv->GetYaxis()->GetTitle());
   THStack *hs = new THStack("hs",title.c_str());

   // Setup the legend
   TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
   l->SetBorderSize(0);
   l->SetNColumns(2);
   l->AddEntry(h_cv,"Central Value","L");

   // Prepare the histograms and add the right information to the legend
   double maximum = 0.0;
   for(size_t i=0;i<hist_v.size();i++){
      if(hist_v.size() > 2){
         hist_v.at(i)->SetLineColor(kGreen);
         if(i == 0) l->AddEntry(hist_v.at(0),"Variations","L");
      }
      else if(hist_v.size() == 1){
         hist_v.at(i)->SetLineColor(kRed);
         l->AddEntry(hist_v.at(0),"Alternative Model","L");
      }
      else if(hist_v.size() == 2 && i == 0){
         hist_v.at(i)->SetLineColor(kRed);
         l->AddEntry(hist_v.at(i),"- 1 #sigma","L");
      }
      else if(hist_v.size() == 2 && i == 1){
         hist_v.at(i)->SetLineColor(kBlue); 
         l->AddEntry(hist_v.at(i),"+ 1 #sigma","L");
      }
      hs->Add(hist_v.at(i));
      if(GetHistMax(hist_v.at(i)) > maximum) maximum = GetHistMax(hist_v.at(i));
   }

   // Prepare the CV histogram
   if(GetHistMax(h_cv) > maximum) maximum = GetHistMax(h_cv);
   h_cv->SetLineColor(1);
   h_cv->SetLineWidth(2);
   h_cv->SetFillColor(0);
   h_cv->SetMaximum(maximum*1.25);
   hs->Add(h_cv);

   // Create the "MicroBooNE" watermark
   TLegend *l_Watermark = new TLegend(0.45,0.900,0.89,0.985);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.005);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetTextFont(62);
   l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");

   // Create the POT label
   TLegend *l_POT = new TLegend(0.54,0.820,0.89,0.900);
   l_POT->SetBorderSize(0);
   l_POT->SetMargin(0.005);
   l_POT->SetTextAlign(32);
   TLegend *l_POT2 = new TLegend(0.54,0.74,0.89,0.82,NULL,"brNDC");
   l_POT2->SetBorderSize(0);
   l_POT2->SetMargin(0.005);
   l_POT2->SetTextAlign(32);

   p_legend->Draw();
   p_legend->cd();
   l->Draw();
   c->cd();
   p_plot->Draw();
   p_plot->cd();

   hs->Draw("HIST nostack");
   hs->GetXaxis()->SetTitleSize(Single_XaxisTitleSize);
   hs->GetYaxis()->SetTitleSize(Single_YaxisTitleSize);
   hs->GetXaxis()->SetTitleOffset(Single_XaxisTitleOffset);
   hs->GetYaxis()->SetTitleOffset(Single_YaxisTitleOffset);
   hs->GetXaxis()->SetLabelSize(Single_XaxisLabelSize);
   hs->GetYaxis()->SetLabelSize(Single_YaxisLabelSize);

   h_cv->Draw("HIST same");

   // Draw the various legends, labels etc.
   if(POT.size() > 0) l_POT->Draw();
   if(POT.size() == 2 && mode.at(0) == kFHC && mode.at(1) == kRHC){
      l_POT2->Draw();
      h_cv->SetMaximum(maximum*1.35);
   }
   if(DrawWatermark) l_Watermark->Draw();

   // If using text labels on the bins
   if(binlabels.size()){
      for(int i=1;i<h_cv->GetNbinsX()+1;i++) hs->GetXaxis()->SetBinLabel(i,binlabels.at(i-1).c_str());
      hs->GetXaxis()->SetLabelSize(Single_TextLabelSize);
   }

   p_plot->Update(); 

   system(("mkdir -p " + plotdir).c_str());
   c->Print((plotdir + "/" + label + "_Sys_" + type + "_" +  dialname + ".png").c_str());
   c->Print((plotdir + "/" + label + "_Sys_" + type + "_" +  dialname + ".pdf").c_str());
   c->Print((plotdir + "/" + label + "_Sys_" + type + "_" +  dialname + ".C").c_str());
   c->Close();

}

///////////////////////////////////////////////////////////////////////////////////////////////

void DrawSystematicBreakdown(TFile* f,TH1D* h_template,vector<string> dials,vector<string> captions,std::string plotdir,std::string label){

   assert(dials.size() == captions.size());

   // Get the covairance matrices and calculate the total covariance of the dials
   std::vector<TH2D*> h_FCov_v(dials.size(),nullptr);
   f->GetObject(("FCov_All_" + dials.at(0)).c_str(),h_FCov_v[0]);
   TH2D* h_FCov = (TH2D*)h_FCov_v.at(0)->Clone("h_FCov");
   for(size_t i_s=1;i_s<dials.size();i_s++){
      f->GetObject(("FCov_All_" + dials.at(i_s)).c_str(),h_FCov_v[i_s]);
      h_FCov->Add(h_FCov_v.at(i_s));      
   }

   // Calculate the fractional error of each dial
   std::vector<TH1D*> h_FracError_v(dials.size(),nullptr);
   for(size_t i_s=0;i_s<dials.size();i_s++){
      h_FracError_v[i_s] = (TH1D*)h_template->Clone(("h_FE_" + dials.at(i_s)).c_str());
      for(int i_b=1;i_b<h_template->GetNbinsX()+1;i_b++) h_FracError_v.at(i_s)->SetBinContent(i_b,sqrt(h_FCov_v.at(i_s)->GetBinContent(i_b,i_b)));
      h_FracError_v.at(i_s)->SetLineColor(GoodLineColors[i_s+1]);
      h_FracError_v.at(i_s)->SetFillStyle(0);
   }

   // Calculate the total fractional error
   TH1D* h_FracError = (TH1D*)h_template->Clone("h_FracError");
   h_FracError->Reset();
  
   double max = 0;  
   for(int i_b=1;i_b<h_template->GetNbinsX()+1;i_b++){
      h_FracError->SetBinContent(i_b,sqrt(h_FCov->GetBinContent(i_b,i_b)));
      max = std::max(max,h_FracError->GetBinContent(i_b));
   }

   // Setup canvas    
   TCanvas* c = new TCanvas("c","c",Single_CanvasX,Single_CanvasY);
   TPad *p_plot = new TPad("p_plot","p_plot",0,0,1,Single_PadSplit);
   TPad *p_legend = new TPad("p_legend","p_legend",0,Single_PadSplit,1,1);
   p_legend->SetBottomMargin(0);
   p_legend->SetTopMargin(0.1);
   p_plot->SetTopMargin(0.01);

   // Create the "MicroBooNE" watermark
   TLegend *l_Watermark = new TLegend(0.45,0.900,0.89,0.985);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.005);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetTextFont(62);
   l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");

   // Draw the beamline plot
   TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
   l->SetBorderSize(0);

   const int nhists = dials.size()+1;
   int ncols = 3;
   if(nhists > 6) ncols = 4;
   if(nhists > 12) ncols = 5;
   l->SetNColumns(ncols);

   // Draw everything
   p_legend->Draw();
   p_legend->cd();
   l->Draw();
   c->cd();
   p_plot->Draw();
   p_plot->cd();

   h_FracError->SetLineColor(GoodLineColors[0]);
   h_FracError->SetLineWidth(3);
   h_FracError->SetLineStyle(2);
   h_FracError->SetFillStyle(0);
   h_FracError->GetYaxis()->SetTitle("Frac. Unc.");
   h_FracError->GetYaxis()->SetRangeUser(0.0,1.15*max);
   h_FracError->Draw("HIST");
   
   for(size_t i_s=0;i_s<dials.size();i_s++){
      h_FracError_v[i_s]->SetLineColor(GoodLineColors[i_s+1]);
      h_FracError_v[i_s]->SetLineWidth(2);
      l->AddEntry(h_FracError_v.at(i_s),captions.at(i_s).c_str(),"L");    
      h_FracError_v[i_s]->Draw("HIST same");
   }

   l->AddEntry(h_FracError,"Total","L");
   h_FracError->Draw("HIST same");

   if(DrawWatermark) l_Watermark->Draw();

   p_plot->RedrawAxis();
      
   system(("mkdir -p " + plotdir).c_str());
   c->Print((plotdir + "/" + label + "_SysBreakdown.pdf").c_str());
   c->Print((plotdir + "/" + label + "_SysBreakdown.C").c_str());
   c->Print((plotdir + "/" + label + "_SysBreakdown.png").c_str());
   c->Close();

   DrawMatrix(h_FCov,h_FCov,label+"_FCov",plotdir,true,true);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void DrawEfficiencyPlot(TEfficiency * Efficiency,std::string title,std::string name,vector<int> mode,vector<double> POT){

   TH1D* h_Before = (TH1D*)Efficiency->GetTotalHistogram()->Clone("h_Before");
   TH1D* h_After = (TH1D*)Efficiency->GetPassedHistogram()->Clone("h_After");

   TCanvas *c = new TCanvas("c1","c1",Single_CanvasX,Single_CanvasY);
   TPad *p_plot = new TPad("pad1","pad1",0,0,1,Single_PadSplit);
   TPad *p_legend = new TPad("pad2","pad2",0,Single_PadSplit,1,1);
   p_legend->SetBottomMargin(0);
   p_legend->SetTopMargin(0.1);
   p_plot->SetTopMargin(0.01);

   TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
   l->SetBorderSize(0);
   l->SetNColumns(2);

   // Create the "MicroBooNE" watermark
   TLegend *l_Watermark = new TLegend(0.45,0.900,0.89,0.985);
   l_Watermark->SetBorderSize(0);
   l_Watermark->SetMargin(0.005);
   l_Watermark->SetTextAlign(32);
   l_Watermark->SetTextSize(0.05);
   l_Watermark->SetTextFont(62);
   l_Watermark->SetHeader("MicroBooNE Simulation, Preliminary","R");

   // Create the POT label
   TLegend *l_POT = new TLegend(0.54,0.820,0.89,0.900);
   l_POT->SetBorderSize(0);
   l_POT->SetMargin(0.005);
   l_POT->SetTextAlign(32);
   TLegend *l_POT2 = new TLegend(0.54,0.74,0.89,0.82,NULL,"brNDC");
   l_POT2->SetBorderSize(0);
   l_POT2->SetMargin(0.005);
   l_POT2->SetTextAlign(32);

   if(mode.at(0) == kFHC) l_POT->SetHeader(("NuMI FHC, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.at(0) == kRHC) l_POT->SetHeader(("NuMI RHC, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.at(0) == kBNB) l_POT->SetHeader(("BNB, " + to_string_with_precision(POT.at(0)/1e20,1) + " #times 10^{20} POT").c_str());

   if(mode.size() == 2 && mode.at(1) == kFHC) l_POT2->SetHeader(("NuMI FHC, " + to_string_with_precision(POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.size() == 2 && mode.at(1) == kRHC) l_POT2->SetHeader(("NuMI RHC, " + to_string_with_precision(POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());
   else if(mode.size() == 2 && mode.at(0) == kBNB && mode.at(1) == kBNB) l_POT->SetHeader(("BNB, " + to_string_with_precision(POT.at(0) + POT.at(1)/1e20,1) + " #times 10^{20} POT").c_str());

   THStack *hs = new THStack("hs",title.c_str());

   h_Before->SetLineWidth(2);
   h_Before->SetLineColor(1);
   hs->Add(h_Before);
   h_After->SetLineWidth(2);
   h_After->SetLineColor(3);
   hs->Add(h_After);

   // Draw everything
   p_legend->Draw();
   p_legend->cd();
   l->Draw();
   c->cd();
   p_plot->Draw();
   p_plot->cd();

   hs->Draw("nostack E0");
   hs->GetXaxis()->SetTitleSize(Single_XaxisTitleSize);
   hs->GetYaxis()->SetTitleSize(Single_YaxisTitleSize);
   hs->GetXaxis()->SetTitleOffset(Single_XaxisTitleOffset);
   hs->GetYaxis()->SetTitleOffset(Single_YaxisTitleOffset);
   hs->GetXaxis()->SetLabelSize(Single_XaxisLabelSize);
   hs->GetYaxis()->SetLabelSize(Single_YaxisLabelSize);

   hs->SetMaximum(1.25*hs->GetMaximum("nostack"));
   if(POT.size() > 1) hs->SetMaximum(hs->GetMaximum("nostack")*1.5);

   p_plot->Update(); 

   // Calculate the efficiencies and their uncertainties
   Efficiency->SetConfidenceLevel(0.68);
   Efficiency->SetStatisticOption(TEfficiency::kBUniform);
   Efficiency->SetPosteriorMode();

   std::vector<double> Efficiency_X;
   std::vector<double> Efficiency_CV;
   std::vector<double> Efficiency_Low;
   std::vector<double> Efficiency_High;

   // Decide on a limit for the RH y axis
   double effmax = 0.0;
   for(size_t i=1;i<h_Before->GetNbinsX()+1;i++)
        effmax = std::max(effmax,Efficiency->GetEfficiency(i)+Efficiency->GetEfficiencyErrorUp(i));

   int binmax = h_Before->GetMaximumBin();
   Float_t rightmax = 1.2*effmax;
   if(POT.size() > 1) rightmax = 1.5*effmax;
   double scale = p_plot->GetUymax()/rightmax;

   for(size_t i=1;i<h_Before->GetNbinsX()+1;i++){
      if(std::isnan(Efficiency->GetEfficiency(i))) continue;
      Efficiency_X.push_back(h_Before->GetBinCenter(i));
      Efficiency_CV.push_back(Efficiency->GetEfficiency(i)*scale);
      Efficiency_Low.push_back(Efficiency->GetEfficiencyErrorLow(i)*scale);
      Efficiency_High.push_back(Efficiency->GetEfficiencyErrorUp(i)*scale);
   }

   TGraphAsymmErrors *g_Efficiency = new  TGraphAsymmErrors(Efficiency_X.size(),&(Efficiency_X[0]),&(Efficiency_CV[0]),0,0,&(Efficiency_Low[0]),&(Efficiency_High[0]));
   g_Efficiency->SetLineColor(kRed);
   g_Efficiency->SetMarkerStyle(5);
   g_Efficiency->SetMarkerSize(2);
   g_Efficiency->SetMarkerColor(kRed);
   g_Efficiency->SetLineWidth(2);

   // Create the new axis object
   TGaxis *axis = new TGaxis(p_plot->GetUxmax(),p_plot->GetUymin(),p_plot->GetUxmax(),p_plot->GetUymax(),0,rightmax,510,"+L");
   axis->SetTitleColor(kRed);
   axis->SetLabelColor(kRed);
   axis->SetTitleSize(Single_YaxisTitleSize);
   axis->SetTitleOffset(0.9*Single_YaxisTitleOffset);
   axis->SetLabelSize(Single_YaxisLabelSize);
   axis->SetTitle("Efficiency");
   axis->Draw();

   l->AddEntry(h_Before,"All Events","L");
   l->AddEntry(g_Efficiency,"Efficiency","P");
   l->AddEntry(h_After,"Selected","L");

   g_Efficiency->Draw("P same");

   // Draw the various legends, labels etc.
   if(POT.size() > 0) l_POT->Draw();
   if(POT.size() == 2 && mode.at(0) == kFHC && mode.at(1) == kRHC) l_POT2->Draw();   
   if(DrawWatermark) l_Watermark->Draw();

   p_plot->Update(); 

   system("mkdir -p Plots/");
   c->Print(("Plots/Efficiency_" + name + "_Ratio.pdf").c_str());
   c->Print(("Plots/Efficiency_" + name + "_Ratio.png").c_str());
   c->Print(("Plots/Efficiency_" + name + "_Ratio.C").c_str());
   c->Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

};

#endif
