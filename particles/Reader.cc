#include "Reader.h"
#include <boost/python.hpp>
using namespace boost::python;
using namespace boost::iostreams;
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>


class ReaderWrap : public Reader, public wrapper<Reader>
{
    public:
        unsigned int getNumParticles() const
        {
            return this->get_override("getNumParticles")();
        }
        unsigned int getNumParticleTypes() const
        {
            return this->get_override("getNumParticleTypes")();
        }
        vec getBox() const
        {
            return this->get_override("getBox")();
        }
        vector< vec >& getRealPos() 
        {
            return this->get_override("getRealPos")();
        }
        vector<vec >& getBoxPos() 
        {
            return this->get_override("getBoxPos")();
        }

};


void export_Reader()
{
    class_<ReaderWrap, boost::shared_ptr<ReaderWrap>, boost::noncopyable>("Reader")
        .def("getNumParticles", pure_virtual(&Reader::getNumParticles))
        .def("getNumParticleTypes", pure_virtual(&Reader::getNumParticleTypes))
        .def("getBox", pure_virtual(&Reader::getBox))
        ;
}
