#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include "../interface/TH1DPlot.hh"


TH1DPlot::~TH1DPlot() {

}


void TH1DPlot::SetMyStyle() {
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

/*void TH1DPlot::CreateCanvas(){


}*/

int TH1DPlot::GeneratePlot() {
    SetMyStyle();
    std::cout << "This is where we would generate the plot" << std::endl;    
    
}
