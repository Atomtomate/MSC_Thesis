#ifndef IOHELPER_H_
#define IOHELPER_H_


#include "GreensFct.hpp"
#include "Config.hpp"

#include "gnuplot-iostream.h"

#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>
#include <string>
#include <memory>
#include <vector>
#include <boost/tuple/tuple.hpp>


#include <boost/filesystem.hpp>

//#include <mgl2/mgl.h>
//#include <mgl2/qt.h>
namespace DMFT
{

struct LogInfos{
    LogInfos(std::string const & _filename): LogInfos(_filename, true, true) {}
    LogInfos(std::string const & _filename, bool printMF, bool printIT):
            filename(_filename), printMF(printMF), printIT(printIT) {}
    std::string filename;
    bool printMF;
    bool printIT;
    std::string infoString = "";
};

class IOhelper//: public mglDraw
{
    public:
        //TODO: thread safe console output
		IOhelper(std::string& _outDir, const Config& c);
        IOhelper(const Config& c);
                
                /*! initialize (new) output dir
                 *  @param [in] _outDir directory path from current working dir
                 */
                void initDir(std::string& _outDir);

		/*! Write all Green's functions to specified folder.
		 *  note that GF must be queued for logging before calling this function.
		 */
		int writeToFile(void) const;
        int writeToFile(const std::string& text, const std::string& file);

        void writeFinalToFile(GreensFct& gf, const LogInfos& li, const bool selfE = false, const RealT U = 0., const bool it_out = false) const;

		/*! Write specific Green's function to file.
		 *  @param  [in]  handle obtained through addGF(const GreensFct* gf, LogInfos l); 
		 */
		int writeToFile(const unsigned index) const;

		/*! Write specific Green's function one time to file.
		 *  @param  [in]  gf    Green's function
         *  @param  [in]  name  filename
		 */
       int writeToFile(GreensFct& gf, std::string & name) const;
                

        void readFromFile(GreensFct& gf, const std::string files_it, const std::string files_mf) const;
        
        /*! Plot Green's function using GnuPlot
         *  @param  [in] gf     Green's function
         */
        static void plot(GreensFct& gf, RealT beta, std::string title);
        static void plot(ImTG& gf, RealT beta, std::string title);
        static void plot(MatG& gf, RealT beta, std::string title);

		/*! Adds Green's function to pool
		 *  @param  gf	constant pointer to Green''s function
		 *  @return	handle for use in writeToFile(const unsigned index)
		 */
		unsigned addGF(GreensFct& gf, LogInfos l) {
			gfList.push_back( std::make_pair(std::ref(gf), l) );
			return gfList.size()-1;
		}

        inline void removeGF(size_t handle)
        {
            gfList.erase(gfList.begin()+handle);
        }


		void setIteration(unsigned i) {iteration = i;}
        //int Draw(mglGraph *gr);
    private:
		using GFListEl = std::pair<std::reference_wrapper<GreensFct>, LogInfos >;
		//TODO: Typename should be some smart pointer (std::shared_ptr), but some performance tests need to be in place first
		std::vector<GFListEl> gfList;
		unsigned iteration;
		boost::filesystem::path outDir;
        const Config& c;

		void writeToFile(GreensFct& gf, const LogInfos& li) const;
};

}
#endif
