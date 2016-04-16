#ifndef MAP2CAT_HPP
#define MAP2CAT_HPP

#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/log/trivial.hpp>
#include <boost/tokenizer.hpp>

#include <healpix_base.h>
#include <healpix_map.h>
#include <pointing.h>
#include <healpix_map_fitsio.h>
#include <datatypes.h>

class Map2Cat
{
public:

    typedef boost::property_tree::ptree propertyTreeType;
    typedef Healpix_Map<double> mapType;
    static const constexpr double deg2rad = M_PI/double(180);
    static const constexpr double rotPhi = double(0);

    explicit Map2Cat(std::string const& iniFileName)
    {
        BOOST_LOG_TRIVIAL(info) << std::string("Reading the ini file ") + std::string(iniFileName);
        boost::property_tree::ini_parser::read_ini(iniFileName,mPropTree);

        std::string inputMapFileName = mPropTree.get<std::string>("input.data_map_file_name");
        BOOST_LOG_TRIVIAL(info) << std::string("Reading the maps from ") + inputMapFileName;
        read_Healpix_map_from_fits(inputMapFileName,mMapN);
        read_Healpix_map_from_fits(inputMapFileName,mMapE1,int(2),int(2));
        read_Healpix_map_from_fits(inputMapFileName,mMapE2,int(3),int(2));
    }

    void generate()
    {
        std::string zBoundsStr = mPropTree.get<std::string>("input.z_bounds");
        BOOST_LOG_TRIVIAL(info) << "z bounds specified as "<< zBoundsStr;
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep(",");
        tokenizer tokens(zBoundsStr, sep);

        for(tokenizer::iterator tokIter = tokens.begin(); tokIter !=tokens.end(); ++tokIter)
        {
            mZBounds.push_back(boost::lexical_cast<double>(*tokIter));
        }

        if(mZBounds.size() != size_t(2))
        {
            std::string msg("The z-bounds should consist of two values. No more, no less.");
            throw std::runtime_error(msg);
        }

        if(mZBounds[0] >= mZBounds[1])
        {
            std::string msg("The upper bound should be greater than the lower bound.");
            throw std::runtime_error(msg);
        }

        size_t numPix = (size_t) mMapN.Npix();

        size_t randomSeed = mPropTree.get<size_t>("input.rand_seed");
        BOOST_LOG_TRIVIAL(info) << "Random seed specified as "<< randomSeed;

        std::mt19937 gen(randomSeed);

        // to generate z values
        std::uniform_real_distribution<double> dist_z(mZBounds[0], mZBounds[1]);
        // to generate e1 and e2 values
        std::normal_distribution<> dist_e(0,1);

        double sigma_e = mPropTree.get<double>("input.sigma_e");
        BOOST_LOG_TRIVIAL(info) << "Std-dvn of ellipticities specified as "<< sigma_e;


        BOOST_LOG_TRIVIAL(info) << "Number of pixel in the map is "<< numPix;

        std::ofstream inputCatFile;
        std::string inputCatFileName = mPropTree.get<std::string>("output.catlogue_file_name");
        BOOST_LOG_TRIVIAL(info) << "Output catalogue file name is "<< inputCatFileName;

        inputCatFile.open( inputCatFileName.c_str(),std::ios::trunc );


        //write the header
        std::string delimiter = mPropTree.get<std::string>("output.delimiter");
        BOOST_LOG_TRIVIAL(info) << "Delimiter for separation is "<< delimiter;

        inputCatFile<<"#"<<"ra"<<delimiter<<"dec"<<delimiter<<"z"
            <<delimiter<<"e1"<<delimiter<<"e2"<<std::endl;


        inputCatFile<<std::setprecision(10);
        // go through the number count map and generate n uniform random from zmin and zmax
        for(size_t i=0;i<numPix;++i)
        {
            size_t numGals = (size_t) mMapN[i];
            double e1 = mMapE1[i];
            double e2 = mMapE2[i];

            pointing pntg = mMapN.pix2ang(i);
            double theta = pntg.theta;
            double phi = pntg.phi;

            double dec = -( theta - M_PI*double(0.5) )/deg2rad;
            double ra = phi/deg2rad + rotPhi;

            for(size_t j=0;j<numGals;++j)
            {
                double zVal = dist_z(gen);
                double e1Val = e1 + dist_e(gen)*sigma_e;
                double e2Val = e2 + dist_e(gen)*sigma_e;
                inputCatFile<<ra<<delimiter<<dec<<delimiter<<zVal<<delimiter<<e1Val<<delimiter<<e2Val<<std::endl;
            }
        }

        inputCatFile.close();
    }

private:
    propertyTreeType mPropTree;
    mapType mMapN;
    mapType mMapE1;
    mapType mMapE2;
    std::vector<double> mZBounds;


};

#endif //MAP2CAT_HPP
