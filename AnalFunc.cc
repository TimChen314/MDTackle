#include "AnalFunc.h"
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/utility.hpp>
using namespace boost::python;

class AnalFuncWrap : public AnalFunc, public wrapper<AnalFunc>
{
    public:
        void compute()  {  this->get_override("compute")(); }
        void finish()  {  this->get_override("finish")(); }
};


void export_AnalFunc()
{
    class_<AnalFuncWrap, boost::shared_ptr<AnalFuncWrap>, boost::noncopyable>("AnalFunc")
        //  .def("compute", pure_virtual(&AnalFunc::compute))
        //  .def("finish", pure_virtual(&AnalFunc::finish))
        ;
}

