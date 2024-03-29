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
        
        THPlot(std::string config_file, std::string postfix="") {
            
            //Fill class member variables from config file
            // Setup options.
            std::string fill_colour;
            po::options_description desc("Options");
            desc.add_options()
              ("name", po::value<std::string>(&name_))
              ("nbinsx", po::value<int>(&nbinsx_))
              ("nbinsy", po::value<int>(&nbinsy_)->default_value(100))
              ("x_start_bin", po::value<double>(&x_start_bin_))
              ("x_end_bin", po::value<double>(&x_end_bin_))
              ("y_start_bin", po::value<double>(&y_start_bin_)->default_value(0.0))
              ("y_end_bin", po::value<double>(&y_end_bin_)->default_value(100.0))
              ("cwidth", po::value<int>(&cwidth_)->default_value(700))
              ("cheight", po::value<int>(&cheight_)->default_value(500))
              ("log_x", po::value<bool>(&log_x_)->default_value(false))
              ("log_y", po::value<bool>(&log_y_)->default_value(false))
              ("output_filename", po::value<std::string>(&output_filename_))
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
            if(twoD_) hist2D_ = TH2D( (name_+postfix).c_str() ,(name_+postfix).c_str() ,nbinsx_,x_start_bin_,x_end_bin_,nbinsy_,y_start_bin_,y_end_bin_);
            else hist_ = TH1D( (name_+postfix).c_str(),(name_+postfix).c_str(),nbinsx_,x_start_bin_,x_end_bin_);
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
        TH1D GetHist() const;
        TH2D Get2DHist() const;
        void SetHist(TH1D const &h);
        void Set2DHist(TH2D const &h);
        std::string GetOutFilename() const;
        void SetOutFilename(std::string const &filename);
        void SetRunInfo(std::pair<std::vector<std::string>,std::vector<std::string> > const &info);
        std::string GetName() const;


    private:
        static void SetMyStyle_();
        double x_start_bin_;
        double x_end_bin_;
        double y_start_bin_;
        double y_end_bin_;
        int nbinsx_;
        int nbinsy_;
        std::string name_;
        bool log_x_;
        bool log_y_;
        bool twoD_;
        std::string x_title_;
        std::string y_title_;
        Int_t r_, g_, b_;
        int cwidth_;
        int cheight_;
        TH1D hist_; 
        TH2D hist2D_; 
        std::string output_filename_;
        std::pair<std::vector<std::string>,std::vector<std::string> > run_info_;

};

#endif
