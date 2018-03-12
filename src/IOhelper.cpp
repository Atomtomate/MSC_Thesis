#include "IOhelper.hpp"

namespace DMFT
{
    IOhelper::IOhelper(std::string& _outDir, const Config& c): c(c)
    {
        initDir(_outDir);
    };
    
    IOhelper::IOhelper(const Config& c): c(c)
    {
        std::string cp = c.outDir;
        initDir(cp);
    };

    void IOhelper::initDir(std::string& _outDir)
    {
        outDir = _outDir;
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
                fs << std::string("solver: \t") << c.solver << std::endl;
                fs << std::string("beta:\t") << c.beta << std::endl;
                fs << std::string("U:   \t") << c.U << std::endl;
                fs << std::string("mu:  \t") << c.mu << std::endl;
                fs << std::string("half BW:  \t") << c.D << std::endl;
                fs << std::string("N_MF: \t") << c.mfCount << std::endl;
                fs << std::string("N_IT: \t") << c.itCount << std::endl;
                fs << std::string("solver: \t") << c.solver << std::endl;
                fs << std::string("additional info: \t") << c.info << std::endl;
                fs.close();
            }
        }
        catch (boost::filesystem::ofstream::failure e)
        {
            LOG(ERROR) << "Error while writing parameters to file: " << e.what();
        }

    }

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
 
    void IOhelper::writeFinalToFile(GreensFct& gf, const LogInfos& li, const bool selfE, const RealT U, const bool it_out) const
    {
        if(!c.isGenerator || c.local.rank() != 0) return ;
        std::string b_str = std::to_string(c.beta);
        std::string u_str = std::to_string(std::max(1.5*U,3.0));
        std::string mu_str = std::to_string(std::min(-1.5*U,-3.0));
        std::string nmf = std::to_string(std::min(_CONFIG_maxMatsFreq, 512));
        std::string file_name = li.filename;
        std::replace(b_str.begin(), b_str.end(), '.', '_');
        //std::replace(u_str.begin(), u_str.end(), '.', '_');
        //std::replace(mu_str.begin(), mu_str.end(), '.', '_');
        std::replace(file_name.begin(), file_name.end(), '.', '_');
        boost::filesystem::path file_maxent = outDir;
        boost::filesystem::path file_maxent_it = outDir;
        boost::filesystem::path file_pade = outDir;
        boost::filesystem::path file_maxent_config = outDir;
        boost::filesystem::path file_pade_config = outDir;
        file_maxent /= boost::filesystem::path( std::string("f_") + file_name + std::string(".out"));
        file_maxent_it /= boost::filesystem::path( std::string("f_") + file_name + std::string(".it.out"));
        file_pade /= boost::filesystem::path( std::string("f_") + file_name + std::string(".pade.out"));
        file_maxent_config /= boost::filesystem::path( std::string("me_") + file_name + std::string("conf.in"));
        file_pade_config /= boost::filesystem::path( std::string("pade_") + file_name + std::string("conf.in"));
        boost::filesystem::ofstream fmt, fmt_it, fmtp, fmt_me_conf, fmt_pade_conf;
        fmt.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        fmt_it.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        fmtp.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        fmt_me_conf.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        fmt_pade_conf.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        const std::string se_add = selfE ? std::string("\nSELF=true\nNORM=")+std::to_string(U*U/16.) : "";
        const std::string me_conf_s = "BETA=" + std::to_string(c.beta) +
                         "\nNDAT=" + nmf +\
                         "\nBACKCONTINUE=false\nMAX_IT=10000" + \
                        "\nNFREQ=1024\nDATASPACE=frequency\nKERNEL=fermionic\nPARTICLE_HOLE_SYMMETRY=true\nDATA=" +
                            std::string("f_") + file_name + std::string(".out") + se_add;
                         //"\nOMEGA_MIN=" + mu_str + "\nOMEGA_MAX=" + u_str +\
        //const std::string pade_conf_s = "imag.BETA=" +std::to_string(c.beta) + "\nimag.NDAT="+ std::to_string(_CONFIG_maxMatsFreq) +\
        //                                 "\nimag.STATISTICS=Fermi\nimag.NEGATIVE_DATA=false\nimag.DATA="+std::string("f_")+li.filename+std::string(".pade.out")+\
        //                                 "\nreal.NFREQ=" + std::to_string(_CONFIG_maxMatsFreq) + "\nreal.OUTPUT="  + li.filename+std::string(".pade.cont")+\
        //                                 "\nreal.FREQUENCY_GRID=linear\npade.PADE_NUMERATOR_DEGREE ="+ std::to_string(_CONFIG_maxMatsFreq/2-1)+\
        //                                 "\npade.PADE_DENOMINATOR_DEGREE =" + std::to_string(_CONFIG_maxMatsFreq/2);
        const std::string pade_conf_s = "imag.BETA=" + std::to_string(c.beta) + "\nimag.NDAT=513\nimag.STATISTICS=Fermi\nimag.NEGATIVE_DATA=false\nimag.DATA="+\
                                         std::string("f_") + file_name + std::string(".pade.out")+\
                                         "\nreal.NFREQ=513\nreal.OUTPUT=" + file_name + std::string(".pade.cont") + \
                                         "\nreal.FREQUENCY_GRID="+ ((U < 4) ? "linear" : "log") +"\nreal.OMEGA_MIN=" + mu_str + "\nreal.OMEGA_MAX=" + u_str +\
                                         "\npade.PADE_NUMERATOR_DEGREE=256\npade.PADE_DENOMINATOR_DEGREE=256";
        try
        {
            fmt.open(file_maxent);
            fmt << gf.getMaxEntString(selfE, 10);
            /*fmt_me_conf.open(file_maxent_config);
            if(it_out)
            {
                fmt_it.open(file_maxent_it);
                fmt_it << gf.getMaxEntItString();
                fmt_it.close();
            }
            fmt_me_conf << me_conf_s;
            fmt.close();
            fmt_me_conf.close();
            //pade:
            fmtp.open(file_pade);
            fmt_pade_conf.open(file_pade_config);
            fmtp << gf.getPadeString();
            fmt_pade_conf << pade_conf_s;
            fmtp.close();
            fmt_pade_conf.close();*/
        }
        catch (boost::filesystem::ofstream::failure e)
        {
            LOG(ERROR) << "Error while writing Green\'s function to file: " << e.what();
        }
    }

    int IOhelper::writeToFile(GreensFct& gf, std::string& name) const
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

    int IOhelper::writeToFile(const std::string& text, const std::string& file)
    {
        if(!c.isGenerator || c.local.rank() != 0) return -1;
        boost::filesystem::path _file = outDir;
        _file	/= boost::filesystem::path(file + std::string(".out"));
        boost::filesystem::ofstream _outStream;
        _outStream.exceptions( std::ofstream::failbit | std::ofstream::badbit );
        try
        {
            _outStream.open(_file);
            _outStream << text;
            _outStream.close();
        }
        catch (boost::filesystem::ofstream::failure e)
        {
            LOG(ERROR) << "Error while writing Matsubara Green\'s function to file: " << e.what();
            LOG(ERROR) << "Error while writing Green\'s function to file: " << e.what();
        }
        return 1;
    }

    void IOhelper::readFromFile(GreensFct& gf, const std::string files_it, const std::string files_mf) const
    {
        boost::filesystem::path file_mf;
        boost::filesystem::path file_it;
        boost::filesystem::ifstream fs;
        fs.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
        if(!c.isGenerator || c.local.rank() != 0) return ;
        if(!files_mf.empty())
        {
            file_mf = boost::filesystem::path(files_mf);
            try
            {
                fs.open(file_mf);
                int n = 0;
                ComplexT up,down;
                fs.close();
            }
            catch (boost::filesystem::ofstream::failure e)
            {
                LOG(ERROR) << "Error while reading Matsubara Green\'s function from file: " << e.what();
            }
        }
        if(!files_it.empty())
        {
            file_it = boost::filesystem::path(files_it);
            try
            {
                fs.open(file_it);
                //TODO: get it oder find highest it
                //TODO: read file into GF
                //TODO: check consistency with fft
                fs.close();
            }
            catch (boost::filesystem::ofstream::failure e)
            {
                LOG(ERROR) << "Error while reading Green\'s function from file: " << e.what();
            }
        }
    }


    void IOhelper::plot(GreensFct& gf, RealT beta, std::string title)
    {
        Gnuplot gp1("gnuplot -persist");
        Gnuplot gp2("gnuplot -persist");
        std::vector<std::pair<RealT,RealT>> mgf_up, mgf_down, itgf_up, itgf_down;
        for(int n=-_CONFIG_maxMatsFreq/2;n<_CONFIG_maxMatsFreq/2;n++)
        {
            mgf_down.push_back( std::make_pair( mFreqS(n, beta), gf.getByMFreq(n,0).imag() ) );
            mgf_up.push_back( std::make_pair( mFreqS(n+_CONFIG_maxMatsFreq/2, beta), gf.getByMFreq(n+_CONFIG_maxMatsFreq/2,1).imag() ) );
        }
        for(RealT i=0;i<beta;i+=beta/_CONFIG_maxTBins)
        {
            itgf_down.push_back( std::make_pair( i, gf(i,0) ) );
            itgf_up.push_back( std::make_pair( i, gf(i,1) ) );
        }    

        gp1 << "set terminal x11\n";
        gp2 << "set terminal x11\n";
        gp1 << "set output ' " << title << "_MF.png'\n";
        gp1 << "plot" << gp1.file1d( mgf_down ) << "with points title ' " <<  title <<  " Matsubara GF, -MF/2 - MF/2', "
            << gp1.file1d( mgf_up ) << "with points title ' " <<  title <<  " Matsubara GF, 0 to MF '" << std::endl;
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
    void IOhelper::plot(MatG& gf, RealT beta, std::string title)
    {
        GreensFct gfTmp(beta);
        gfTmp.setByMFreq(gf);
        gfTmp.transformMtoT();
        plot(gfTmp, beta, title);
    }

}
/*int IOhelper::Draw(mglGraph *gr)
  {
  return 0;
  }*/
