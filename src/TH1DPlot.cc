#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TColor.h>
#include <iostream>
#include "../interface/TH1DPlot.hh"


//Private functions

void TH1DPlot::SetMyStyle_() {
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

    // color palette for 2D temperature plots
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


std::string TH1DPlot::GetName(){
    return name;
}

TH1D* TH1DPlot::GetHist(){
    return hist;
}

void TH1DPlot::SetHist(TH1D* h){
    hist = h;
}

std::string TH1DPlot::GetOutFilename(){
    return output_filename;
}

void TH1DPlot::SetOutFilename(std::string filename){
    output_filename = filename;
}

void TH1DPlot::SetRunInfo(std::pair<std::string,std::string> info){
    run_info = info;
}


int TH1DPlot::GeneratePlot() {
    
    SetMyStyle_();
    TCanvas* c1 = new TCanvas();
    hist->SetFillStyle(1001);
    hist->SetFillColor(TColor::GetColor(r_,g_,b_));
    hist->GetYaxis()->SetTitle( y_title_.c_str() );
    hist->GetXaxis()->SetTitle( x_title_.c_str() );
    hist->SetMinimum(0.0);
    hist->Draw("hist");
    TLatex *title_latex = new TLatex();
    title_latex->SetNDC();
    title_latex->SetTextSize(0.04);
    title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str() );
    c1->SaveAs((output_filename+".png").c_str());
    c1->SaveAs((output_filename+".pdf").c_str());
    delete c1;
    
    return 0; 
}
