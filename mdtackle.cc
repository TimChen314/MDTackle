#include "XmlReader.h"
#include "Traj2CMBin.h"
#include "BinReader.h"
#include "XmlModify.h"
#include "BasicInfo.h"
#include "MolInfo.h"
#include "Traj.h"
#include "AnalFunc.h"
#include "SegRgDis.h"
#include "Analysis.h"
#include "BoxLen.h"
#include "RgByType.h"
#include<boost/python.hpp>
using namespace boost::python;

void OutputVersionInfo()
{
    cout<<"Version test 1.0 "<<endl;
}


BOOST_PYTHON_MODULE(mdtackle) // It's important that the library file is named like you declare the module here: BOOST_PYTHON_MODULE(hello_ext) that is hello_ext.dll or hello_ext.so
{
    OutputVersionInfo();
    export_BasicInfo();
    export_Traj();
    export_Reader();
    export_XmlReader();
    export_BinReader();
    export_XmlModify();
    export_AnalFunc();
    export_SegRgDis();
    export_Analysis();
    export_BoxLen();
    export_Traj2CMBin();
    export_RgByType();

}


