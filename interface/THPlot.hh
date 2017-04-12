#ifndef PlotProcessedData_THPlot_hh
#define PlotProcessedData_THPlot_hh

#include <vector>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/program_options.hpp>


namespace po = boost::program_options;

class THPlot {


    public:

        THPlot(){;}
        
        THPlot(std::string config_file) {
            
            //Fill class member variables from config file
            // Setup options.
            std::string fill_colour;
            po::options_description desc("Options");
            desc.add_options()
              ("name", po::value<std::string>(&name))
              ("nbinsx", po::value<int>(&nbinsx))
              ("nbinsy", po::value<int>(&nbinsy)->default_value(100))
              ("x_start_bin", po::value<double>(&x_start_bin_))
              ("x_end_bin", po::value<double>(&x_end_bin_))
              ("y_start_bin", po::value<double>(&y_start_bin_)->default_value(0.0))
              ("y_end_bin", po::value<double>(&y_end_bin_)->default_value(100.0))
              ("cwidth", po::value<int>(&cwidth)->default_value(700))
              ("cheight", po::value<int>(&cheight)->default_value(500))
              ("log_x", po::value<bool>(&log_x_)->default_value(false))
              ("log_y", po::value<bool>(&log_y_)->default_value(false))
              ("output_filename", po::value<std::string>(&output_filename))
              ("x_title", po::value<std::string>(&x_title_))
              ("y_title", po::value<std::string>(&y_title_))
              ("twoD", po::value<bool>(&twoD_)->default_value(false))
              ("fill_colour", po::value<std::string>(&fill_colour)->default_value("0,0,0"));
            
            po::variables_map vm;
            po::notify(vm);

            // Write, read, and print settings.
            std::ifstream settings_file((config_file).c_str());

            // Clear the map.
            vm = po::variables_map();

            po::store(po::parse_config_file<char>(settings_file , desc), vm);
            po::notify(vm);  
            
            SetMyStyle_();
            if(twoD_) hist2D = new TH2D(name.c_str(),name.c_str(),nbinsx,x_start_bin_,x_end_bin_,nbinsy,y_start_bin_,y_end_bin_);
            else hist = new TH1D(name.c_str(),name.c_str(),nbinsx,x_start_bin_,x_end_bin_);
            std::stringstream ss(fill_colour);
            std::vector<Int_t> vect;
            Float_t i;
            while (ss >> i)
            {
              vect.push_back(i);
              if (ss.peek() == ',')
               ss.ignore();
            }
            r_ = vect[0];
            g_ = vect[1];
            b_ = vect[2];
        }
        
        ~THPlot(){
        };
        
        
        int GeneratePlot();
        TH1D* GetHist();
        TH2D* Get2DHist();
        void SetHist(TH1D* h);
        void Set2DHist(TH2D* h);
        std::string GetOutFilename();
        void SetOutFilename(std::string filename);
        void SetRunInfo(std::pair<std::string,std::string> info);

        std::string name;
        int nbinsx;
        int nbinsy;
        int cwidth;
        int cheight;
        TH1D* hist; 
        TH2D* hist2D; 
        std::string output_filename;
        std::pair<std::string,std::string> run_info;
        std::string GetName();

    private:
        int *ptr;
        static void SetMyStyle_();
        double x_start_bin_;
        double x_end_bin_;
        double y_start_bin_;
        double y_end_bin_;
        bool log_x_;
        bool log_y_;
        bool twoD_;
        std::string x_title_;
        std::string y_title_;
        Int_t r_, g_, b_;

};

#endif
