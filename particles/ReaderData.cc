#include <ReaderData.h>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace boost;
using namespace boost::iostreams;
using namespace boost::python;


ReaderData::ReaderData(vec& para_box,vector<vec>& para_pos,vector<vec_int>& para_image, vector<vec>& para_vel, \
        vector<unsigned>& para_type,vector<double>& para_mass, vector<Bond>& para_bond, vector<Angle>& para_angle, \
        vector < string > & para_type_mapping)
{
    m_box=para_box;
    m_pos_array=para_pos;
    m_image_array=para_image;
    m_vel_array=para_vel;
    m_idx_type_array=para_type;
    m_type_mapping = para_type_mapping;
    m_mass_array=para_mass;
    m_bonds=para_bond;
    m_angles=para_angle;

    if(m_type_mapping.size()!=0 && m_idx_type_array.size()!=0)
        for(auto v : m_idx_type_array)
            m_string_type_array.push_back(m_type_mapping[v]);
    else{
        cerr<<"*****ReaderData: Can't get type array (the type of type_array is string). *****"<<endl;
        throw runtime_error("ReaderData: Can't get type array (the type of type_array is string).");
    }

    m_ExistTerm.insert("box");
    if(m_pos_array.size()!=0) m_ExistTerm.insert("position");
    if(m_image_array.size()!=0) m_ExistTerm.insert("image");
    if(m_vel_array.size()!=0) m_ExistTerm.insert("velocity");
    if(m_string_type_array.size()!=0) m_ExistTerm.insert("type");
    if(m_mass_array.size()!=0) m_ExistTerm.insert("mass");
    if(m_bonds.size()!=0) m_ExistTerm.insert("bond");
    if(m_angles.size()!=0) m_ExistTerm.insert("angle");

}


vec ReaderData::getBox() const
{
    return m_box;
}
unsigned int ReaderData::getNumParticles() const
{
    return m_pos_array.size();
}
unsigned int ReaderData::getNumParticleTypes() const
{
    return (unsigned int)m_type_mapping.size();
}
vector< vec >& ReaderData::getRealPos() 
{
    m_realpos_array=m_pos_array;
    return m_realpos_array;
}
vector<vec >& ReaderData::getBoxPos() 
{
    return m_pos_array;
}





