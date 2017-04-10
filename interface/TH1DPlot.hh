#ifndef PlotProcessedData_TH1DPlot_h
#define PlotProcessedData_TH1DPlot_h

#include <vector>
#include <string>
#include <TH1D.h>



class TH1DPlot {


    public:

        TH1DPlot(){
            nbins = 100;
            x_start_bin = 0;
            x_end_bin = 100;
            log_x=false;
            log_y=false;
        }
        ~TH1DPlot();
        int GeneratePlot();

        int nbins;
        double x_start_bin;
        double x_end_bin;
        bool log_x;
        bool log_y;
        std::string output_filename;
        std::string title_left;
        static void SetMyStyle();
        


    private:


};

#endif
