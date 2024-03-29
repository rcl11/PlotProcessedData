#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TColor.h>
#include <iostream>
#include <string>
#include "../interface/THPlot.hh"


//Private functions

void THPlot::SetMyStyle_() {
    TStyle *MyStyle = new TStyle("My-Style","Stealing some cosmetics from Htt!");
    gStyle = MyStyle;

    // Canvas
    MyStyle->SetCanvasColor     (0);
    MyStyle->SetCanvasBorderSize(10);
    MyStyle->SetCanvasBorderMode(0);
    MyStyle->SetCanvasDefH      (600);
    MyStyle->SetCanvasDefW      (600);
    MyStyle->SetCanvasDefX      (100);
    MyStyle->SetCanvasDefY      (100);

    // color palette for twoD temperature plots
    MyStyle->SetPalette(1,0);
    
    // Pads
    MyStyle->SetPadColor       (0);
    MyStyle->SetPadBorderSize  (10);
    MyStyle->SetPadBorderMode  (0);
    MyStyle->SetPadBottomMargin(0.13);
    MyStyle->SetPadTopMargin   (0.08);
    MyStyle->SetPadLeftMargin  (0.17);
    MyStyle->SetPadRightMargin (0.05);
    MyStyle->SetPadGridX       (0);
    MyStyle->SetPadGridY       (0);
    MyStyle->SetPadTickX       (1);
    MyStyle->SetPadTickY       (1);

    // Frames
    MyStyle->SetLineWidth(2);
    MyStyle->SetFrameFillStyle ( 0);
    MyStyle->SetFrameFillColor ( 0);
    MyStyle->SetFrameLineColor ( 1);
    MyStyle->SetFrameLineStyle ( 0);
    MyStyle->SetFrameLineWidth ( 2);
    MyStyle->SetFrameBorderSize(10);
    MyStyle->SetFrameBorderMode( 0);

    // Histograms
    MyStyle->SetHistFillColor(2);
    MyStyle->SetHistFillStyle(0);
    MyStyle->SetHistLineColor(1);
    MyStyle->SetHistLineStyle(0);
    MyStyle->SetHistLineWidth(3);
    MyStyle->SetNdivisions(510);

    // Functions
    MyStyle->SetFuncColor(1);
    MyStyle->SetFuncStyle(0);
    MyStyle->SetFuncWidth(2);

    // Various
    MyStyle->SetMarkerStyle(20);
    MyStyle->SetMarkerColor(kBlack);
    MyStyle->SetMarkerSize (1.1);

    MyStyle->SetTitleBorderSize(0);
    MyStyle->SetTitleFillColor (0);
    MyStyle->SetTitleX         (0.2);

    MyStyle->SetTitleSize  (0.055,"X");
    MyStyle->SetTitleOffset(1.200,"X");
    MyStyle->SetLabelOffset(0.005,"X");
    MyStyle->SetLabelSize  (0.040,"X");
    MyStyle->SetLabelFont  (42   ,"X");

    MyStyle->SetStripDecimals(kFALSE);
    MyStyle->SetLineStyleString(11,"20 10");
    MyStyle->SetTitleSize  (0.055,"Y");
    MyStyle->SetTitleOffset(1.600,"Y");
    MyStyle->SetLabelOffset(0.010,"Y");
    MyStyle->SetLabelSize  (0.040,"Y");
    MyStyle->SetLabelFont  (42   ,"Y");

    MyStyle->SetOptTitle(0);

    MyStyle->SetTextSize   (0.055);
    MyStyle->SetTextFont   (42);

    MyStyle->SetStatFont   (42);
    MyStyle->SetTitleFont  (42);
    MyStyle->SetTitleFont  (42,"X");
    MyStyle->SetTitleFont  (42,"Y");
    MyStyle->SetEndErrorSize(0);

    MyStyle->SetOptStat    (0);

    gROOT->ForceStyle();
    return;


}

//Public functions


std::string THPlot::GetName() const{
    return name_;
}

TH1D THPlot::GetHist() const {
    return hist_;
}

TH2D THPlot::Get2DHist() const {
    return hist2D_;
}

void THPlot::SetHist(TH1D const &h){
    hist_ = h;
}

void THPlot::Set2DHist(TH2D const &h){
    hist2D_ = h;
}

std::string THPlot::GetOutFilename() const {
    return output_filename_;
}

void THPlot::SetOutFilename(std::string const &filename){
    output_filename_ = filename;
}

void THPlot::SetRunInfo(std::pair<std::vector<std::string>,std::vector<std::string> > const &info){
    run_info_ = info;
}


int THPlot::GeneratePlot() {
    
    SetMyStyle_();
    TCanvas* c1 = new TCanvas("c1","c1",cwidth_,cheight_);
    if(!twoD_){
      if(hist_.GetSumw2N() == 0) hist_.Sumw2();
      hist_.GetYaxis()->SetTitle( y_title_.c_str() );
      hist_.GetXaxis()->SetTitle( x_title_.c_str() );
      hist_.SetFillStyle(1001);
      hist_.SetFillColor(TColor::GetColor(r_,g_,b_));
      hist_.SetMinimum(0.0);
      hist_.Draw("hist");
    } else {
      hist2D_.GetYaxis()->SetTitle( y_title_.c_str() );
      hist2D_.GetXaxis()->SetTitle( x_title_.c_str() );
      hist2D_.SetMarkerStyle(20);
      hist2D_.Draw("colz");
    }
    TLatex *title_latex = new TLatex();
    title_latex->SetNDC();
    title_latex->SetTextSize(0.04);
    std::string runs="";
    std::string subruns="";
    if(run_info_.first.size() > 1){
      subruns = "all";  
      for(unsigned i=0; i<run_info_.first.size(); i++){
        if(i==0) {
            //remove preceding r and 0s
            runs.append(((run_info_.first)[i]).erase(0,(run_info_.first)[i].find_first_not_of("r0")));
            runs.append("-");
        }
        if(i==run_info_.first.size()-1) {
            runs.append(((run_info_.first)[i]).erase(0,(run_info_.first)[i].find_first_not_of("r0")));
        }
      }
    } else {
      runs = ((run_info_.first)[0]).erase(0,(run_info_.first)[0].find_first_not_of("r0"));  
      for(unsigned i=0; i<run_info_.second.size(); i++){
         //remove preceding s and 0s
         subruns.append(((run_info_.second)[i]).erase(0,2));
         if(i!=run_info_.second.size()-1) subruns.append(",");
      }
    }
    double start_pos =0.5;
    if(twoD_) start_pos=0.35;
    title_latex->DrawLatex(start_pos, 0.94, ("Runs: "+ runs + " Subruns: " + subruns).c_str() );
    if(twoD_) c1->SetRightMargin(0.15);
    c1->Update();
    c1->SaveAs((output_filename_+".png").c_str());
    c1->SaveAs((output_filename_+".pdf").c_str());
    delete c1;
    delete title_latex;
    
    return 0; 
}
