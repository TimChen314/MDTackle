#include "MolDataStruct.h"
#include "Reader.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <time.h>

#include <iomanip>
using namespace std;

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>


#ifndef __DATA_BULDER_H__
#define __DATA_BULDER_H__

class ReaderData : public Reader
{
    public:
        ReaderData(vec& para_box,vector<vec>& para_pos,vector<vec_int>& para_image, vector<vec>& para_vel, \
                vector<unsigned>& para_type,vector<double>& para_mass, vector<Bond>& para_bond,   \
                vector<Angle>& para_angle, vector < string > & para_type_mapping); 
        // NOTE!!! : if you don't want to copy the data, you can code a version which get 
        //           data by swap(). In this way, old data will be destoryed.
        virtual vec getBox() const;
        virtual unsigned int getNumParticles() const;
        virtual unsigned int getNumParticleTypes() const;
        virtual vector< vec >& getRealPos() ;
        virtual vector<vec >& getBoxPos() ;

    private:

};


#endif
