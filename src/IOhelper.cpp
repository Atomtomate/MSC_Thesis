#include "IOhelper.hpp"

namespace DMFT
{
    IOhelper::IOhelper(const std::string& _outDir, const Config& c): outDir(_outDir), c(c){
        if(!c.isGenerator || !c.local.rank() == 0) return;
        try
        {
            if(!boost::filesystem::create_directory(outDir))
            {
                if(!boost::filesystem::is_directory(outDir))
                    LOG(ERROR) << "Output directory for IOhelper could not be created!";
            }


        }
        catch(boost::filesystem::filesystem_error const & err)
        {
            LOG(ERROR) << err.what();
        }
        try
        {
            if(c.isGenerator && c.local.rank() == 0)
            {
                boost::filesystem::path file = outDir;
                file /= boost::filesystem::path(std::string("params.out"));
                boost::filesystem::ofstream fs;
                fs.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
                fs.open (file);
                fs << std::string("beta:\t") << c.beta << std::endl;
                fs << std::string("U:   \t") << c.U << std::endl;
                fs << std::string("mu:  \t") << c.mu << std::endl;
                fs.close();
            }
        }
        catch (boost::filesystem::ofstream::failure e)
        {
            LOG(ERROR) << "Error while writing parameters to file: " << e.what();
        }
    };

    int IOhelper::writeToFile(void) const
    {
        if(!c.isGenerator || c.local.rank() != 0) return 0;
        for(auto el: gfList) writeToFile(std::get<0>(el).get(), std::get<1>(el));
        return 0;
    }

    int IOhelper::writeToFile(unsigned index) const
    {
        if(!c.isGenerator || c.local.rank() != 0) return 0;
        if(index >= gfList.size()) return 1;	// error, index out of range
        writeToFile(std::get<0>(gfList[index]), std::get<1>(gfList[index]));
        return 0;
    }

    int IOhelper::writeToFile(GreensFct& gf, std::string & name) const
    {
        if(!c.isGenerator || c.local.rank() != 0) return 0;
        LogInfos tmp(name);
        writeToFile(gf, tmp);
        return 0;
    }

    void IOhelper::writeToFile(GreensFct& gf,const LogInfos& li) const
    {
        if(!c.isGenerator || c.local.rank() != 0) return ;
        boost::filesystem::path file_mf = outDir;
        boost::filesystem::path file_it = outDir;
        boost::filesystem::path file_maxent = outDir;
        file_mf	/= boost::filesystem::path(std::to_string(iteration) + std::string("_") + li.filename + std::string("_MF.out"));
        file_it	/= boost::filesystem::path(std::to_string(iteration) + std::string("_") + li.filename + std::string("_IT.out"));
        file_maxent /= boost::filesystem::path(std::to_string(iteration) + std::string("_") + li.filename + std::string("_MaxEnt.out"));
        boost::filesystem::ofstream fm,ft,fmt;
        fm.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        ft.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        fmt.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        try
        {
            fm.open (file_mf);
            ft.open (file_it);
            fmt.open(file_maxent);
            fm << li.infoString;
            ft << li.infoString;
            fm << gf.getMGFstring();
            ft << gf.getITGFstring();
            fmt << gf.getMaxEntString();
            fm.close(); ft.close();fmt.close();
        }
        catch (boost::filesystem::ofstream::failure e)
        {
            LOG(ERROR) << "Error while writing Matsubara Green\'s function to file: " << e.what();
            LOG(ERROR) << "Error while writing Green\'s function to file: " << e.what();
        }
    }


    void IOhelper::plot(GreensFct& gf, RealT beta, std::string title)
    {
        Gnuplot gp1("gnuplot -persist");
        Gnuplot gp2("gnuplot -persist");
        std::vector<std::pair<RealT,RealT>> mgf_up, mgf_down, itgf_up, itgf_down;
        for(int n=0;n<_CONFIG_maxMatsFreq;n++)
        {
            mgf_down.push_back( std::make_pair( mFreq(n, beta), gf.getByMFreq(n,0).imag() ) );
            mgf_up.push_back( std::make_pair( mFreq(n, beta), gf.getByMFreq(n,1).imag() ) );
        }
        for(RealT i=0;i<beta;i+=beta/_CONFIG_maxTBins)
        {
            itgf_down.push_back( std::make_pair( i, gf(i,0) ) );
            itgf_up.push_back( std::make_pair( i, gf(i,1) ) );
        }    

        gp1 << "set output ' " << title << "_MF.png'\n";
        gp1 << "plot" << gp1.file1d( mgf_down ) << "with points title ' " <<  title <<  " Matsubara GF, sigma down', "
            << gp1.file1d( mgf_up ) << "with points title ' " <<  title <<  " Matsubara GF, sigma up'" << std::endl;
        gp2 << "set output ' " << title << "_IT.png'\n";
        gp2 << "plot" << gp2.file1d( itgf_up) << "with points title ' " <<  title <<  " iTime GF, sigma down', "
            << gp2.file1d( itgf_down ) << "with points title ' " <<  title <<  " iTime GF. sigma up'" << std::endl;
    }

    void IOhelper::plot(ImTG& gf, RealT beta, std::string title)
    {
        GreensFct gfTmp(beta);
        gfTmp.setByT(gf);
        gfTmp.transformTtoM();
        plot(gfTmp, beta, title);
    }

}
/*int IOhelper::Draw(mglGraph *gr)
  {
  return 0;
  }*/
