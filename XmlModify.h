#include "XmlReader.h"
#include "ReaderData.h"
#include "DumpXml.h"
#include "BasicInfo.h"
#include "TopoInfo.h"
#include <set>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <boost/python.hpp>
using namespace std;

#ifndef __MODIFY_h__
#define __MODIFY_h__

class XmlModify
{
    public: 
        XmlModify(boost::shared_ptr<Reader> read);

        boost::shared_ptr<Reader> build;

        void XmlCutTail();
        void XmlCutOutOneSphere(unsigned na_sphere);
        void XmlMagnifyBox(BasicInfo&, unsigned, unsigned, unsigned);
        //  void XmlMagnifyBox();
    private:
        vec box;
        vector<vec> pos;
        vector<vec> vel;
        vector<vec_int> image;
        vector<unsigned int> type;
        vector<double> mass;
        vector<Bond> bond;
        vector<Angle> angle;
        vector<string> typemapping;
        set<string> ExistTerm;
};

void export_XmlModify();

#endif
