#include "BinReader.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace boost::iostreams;
using namespace boost::python;

BinReader::BinReader(const std::string& fname)
{
    // initialize member variables
    m_timestep = 0;
    m_num_dimensions = 3;
    file_version = -1;
    // read in the file
    readFile(fname);
    convertToVec();
}

BinReader::BinReader(const std::string& fname, bool m_only_parse_head_tmp)
    : m_only_parse_head(m_only_parse_head_tmp)
{
    m_timestep = 0;
    m_num_dimensions = 3;
    file_version = -1;
    readFile(fname);
}

BinReader::BinReader(const std::string& fname, const unsigned block_begin_, const unsigned block_end_)
    : block_begin(block_begin_)
    , block_end(block_end_)
{
    m_timestep = 0;
    m_num_dimensions = 3;
    file_version = -1;
    readFile(fname);
}
void BinReader::ReadExistTerm()
{
    m_ExistTerm.insert("box");
    if (m_input_position)
        m_ExistTerm.insert("position");
    if (m_input_velocity)
        m_ExistTerm.insert("velocity");
    if (m_input_image)
        m_ExistTerm.insert("image");
    if (m_input_type)
        m_ExistTerm.insert("type");
    if (m_input_mass)
        m_ExistTerm.insert("mass");
    if (m_input_diameter)
        m_ExistTerm.insert("diameter");
    if (m_input_body)
        m_ExistTerm.insert("body");
    if (m_input_accel)
        m_ExistTerm.insert("accel");
    if (m_input_charge)
        m_ExistTerm.insert("charge");
    if (m_input_moleculeid)
        m_ExistTerm.insert("moleculeid");
    if (m_input_virial)
        m_ExistTerm.insert("virial");
    if (m_input_force)
        m_ExistTerm.insert("force");
    if (m_input_anisotropy)
        m_ExistTerm.insert("anisotropy");
    if (m_input_bond)
        m_ExistTerm.insert("bond");
    if (m_input_angle)
        m_ExistTerm.insert("angle");
    if (m_input_dihedral)
        m_ExistTerm.insert("dihedral");
    if (m_input_integrator)
        m_ExistTerm.insert("integrator");
    if (m_input_rigid)
        m_ExistTerm.insert("rigid");
    if (m_input_cris)
        m_ExistTerm.insert("cris");
    if (m_input_init)
        m_ExistTerm.insert("init");
}

void BinReader::convertToVec()
{
    if (file_version == 3201) {
        for (unsigned i = 0; i != m_x_array.size(); ++i) {
            vec tmpv(m_x_array[i], m_y_array[i], m_z_array[i]);
            m_pos_array.push_back(tmpv);
        }
        for (unsigned i = 0; i != m_ix_array.size(); ++i) {
            vec_int tmpv(m_ix_array[i], m_iy_array[i], m_iz_array[i]);
            m_image_array.push_back(tmpv);
        }
        for (unsigned i = 0; i != m_vx_array.size(); ++i) {
            vec tmpv(m_vx_array[i], m_vy_array[i], m_vz_array[i]);
            m_vel_array.push_back(tmpv);
        }
    } else if (file_version == 3202) {
        for (unsigned i = 0; i != m_x_array.size(); ++i) {
            vec tmpv(m_x_array[i], m_y_array[i], m_z_array[i]);
            m_realpos_array.push_back(tmpv);
        }
        for (unsigned i = 0; i != m_vx_array.size(); ++i) {
            vec tmpv(m_vx_array[i], m_vy_array[i], m_vz_array[i]);
            m_vel_array.push_back(tmpv);
        }
    } else if (file_version == 1) {
        //has tackled in readFile()
        if (m_input_type || m_input_charge || m_input_moleculeid || m_input_virial || m_input_force
            || m_input_anisotropy || m_input_init || m_input_cris || m_input_bond || m_input_angle
            || m_input_dihedral || m_input_integrator)
            cout << "Error! This program can't tackle info except pos, image and vec." << endl;
    } else if (file_version == 2) {
        //has tackled in readFile()
        if (m_input_type || m_input_charge || m_input_body || m_input_virial || m_input_force
            || m_input_anisotropy || m_input_init || m_input_cris || m_input_bond || m_input_angle
            || m_input_dihedral || m_input_integrator)
            cout << "Error! This program can't tackle info except pos, image and vec." << endl;
    } else if (file_version == 3) {
        //has tackled in readFile()
        if (m_input_type || m_input_mass || m_input_diameter || m_input_body || m_input_accel || m_input_charge || m_input_moleculeid || m_input_virial || m_input_force || m_input_anisotropy || m_input_bond || m_input_angle || m_input_dihedral || m_input_integrator || m_input_rigid)
            cout << "Error! This program can't tackle info except pos, image and vec." << endl;
    } else {
        cerr << "*****Error! file_version does not exist.*****" << endl;
        throw runtime_error("Error! file_version does not exist.");
    }
}

vector<vec>& BinReader::getRealPos()
{
    if (file_version == 3201 || file_version == 1 || file_version == 2 || file_version == 3) {
        for (unsigned i = 0; i != m_pos_array.size(); ++i)
            m_realpos_array.push_back(m_pos_array[i] += m_box * m_image_array[i]);
    } else if (file_version != 3202) { //if file_version == 3202, then nothing should be done.
        cerr << "***** Error! Unknown file type *****" << endl;
        throw runtime_error("Error! Unknown file type ");
    }
    return m_realpos_array;
}

vector<vec>& BinReader::getBoxPos()
{
    if (file_version != 3202) {
        vec boxINV(vec(1.0, 1.0, 1.0) / m_box);
        m_pos_array = m_realpos_array;
        for (unsigned i = 0; i != m_realpos_array.size(); ++i)
            pbc(m_box, boxINV, m_pos_array[i]);
    } else if (!(file_version == 3201 || file_version == 1 || file_version == 2 || file_version == 3)) {
        cerr << "***** Error! Unknown file type *****" << endl;
        throw runtime_error("Error! Unknown file type ");
    }
    return m_pos_array;
}

unsigned int BinReader::getNumDimensions() const
{
    return m_num_dimensions;
}

unsigned int BinReader::getNumParticles() const
{
    assert(m_x_array.size() > 0);
    return (unsigned int)m_x_array.size();
}

unsigned int BinReader::getNumParticleTypes() const
{
    assert(m_type_mapping.size() > 0);
    return (unsigned int)m_type_mapping.size();
}

vec BinReader::getBox() const
{
    return m_box;
}

unsigned int BinReader::getTimeStep() const
{
    return m_timestep;
}

/* change internal timestep number. */
void BinReader::setTimeStep(unsigned int ts)
{
    m_timestep = ts;
}

//! Helper function to read a string from the file
static string read_string(istream& f)
{
    unsigned int len;
    f.read((char*)&len, sizeof(unsigned int));
    if (len != 0) {
        char* cstr = new char[len + 1];
        f.read(cstr, len * sizeof(char));
        cstr[len] = '\0';
        string str(cstr);
        delete[] cstr;
        return str;
    } else
        return string();
}

void BinReader::readFile(const string& fname)
{
    // check to see if the file has a .gz extension or not and enable decompression if it is
    bool enable_decompression = false;
    string ext = fname.substr(fname.size() - 3, fname.size());
    if (ext == string(".gz"))
        enable_decompression = true;

    filtering_istream f;
    if (enable_decompression)
        f.push(gzip_decompressor());
    f.push(file_source(fname.c_str(), ios::in | ios::binary));

    // handle errors
    if (f.fail()) {
        cerr << endl
             << "***Error! Error opening " << fname << endl
             << endl;
        throw runtime_error("Error reading binary file");
    }

    // read magic
    unsigned int magic = 0x444d4f48;
    unsigned int file_magic;
    f.read((char*)&file_magic, sizeof(int));
    if (magic != file_magic) {
        cerr << endl
             << "***Error! " << fname << " does not appear to be a _bin file." << endl;
        if (enable_decompression)
            cerr << "Is it perhaps an uncompressed file with an erroneous .gz extension?" << endl
                 << endl;
        else
            cerr << "Is it perhaps a compressed file without a .gz extension?" << endl
                 << endl;

        throw runtime_error("Error reading binary file");
    }

    f.read((char*)&file_version, sizeof(int));
    // right now, the version tag doesn't do anything: just warn if they don't match
    if (file_version != 3202 && m_read_by_block) {
        cerr << "Error! Read by block can only be used while the file version is 3202." << endl
             << "Version of present file " << fname << " is " << file_version << ". " << endl;
        throw runtime_error("Error Read by block");
    }
    if ((enable_decompression) && m_read_by_block) {
        cerr << "Error! Read by block can only be used while the file is not gziped." << endl;
        throw runtime_error("Error Can not Read by block");
    }

    if (file_version == 3201) {
        f.read((char*)&m_input_position, sizeof(bool));
        f.read((char*)&m_input_image, sizeof(bool));
        f.read((char*)&m_input_velocity, sizeof(bool));
        //ct: add at 2014.03.08
        ReadExistTerm();

        int timestep;
        f.read((char*)&timestep, sizeof(unsigned int));
        m_timestep = timestep;
        unsigned int dimensions;
        f.read((char*)&dimensions, sizeof(unsigned int));
        m_num_dimensions = dimensions;
        float Lx, Ly, Lz;
        f.read((char*)&Lx, sizeof(float));
        f.read((char*)&Ly, sizeof(float));
        f.read((char*)&Lz, sizeof(float));
        m_box = vec(Lx, Ly, Lz);
        unsigned int np = 0;
        f.read((char*)&np, sizeof(unsigned int));
        if (m_only_parse_head)
            return;

        if (m_input_position) {
            m_x_array.resize(np);
            m_y_array.resize(np);
            m_z_array.resize(np);
            f.read((char*)&(m_x_array[0]), np * sizeof(float));
            f.read((char*)&(m_y_array[0]), np * sizeof(float));
            f.read((char*)&(m_z_array[0]), np * sizeof(float));
        }
        if (m_input_image) {
            m_ix_array.resize(np);
            m_iy_array.resize(np);
            m_iz_array.resize(np);
            f.read((char*)&(m_ix_array[0]), np * sizeof(int));
            f.read((char*)&(m_iy_array[0]), np * sizeof(int));
            f.read((char*)&(m_iz_array[0]), np * sizeof(int));
        }
        if (m_input_velocity) {
            m_vx_array.resize(np);
            m_vy_array.resize(np);
            m_vz_array.resize(np);
            f.read((char*)&(m_vx_array[0]), np * sizeof(float));
            f.read((char*)&(m_vy_array[0]), np * sizeof(float));
            f.read((char*)&(m_vz_array[0]), np * sizeof(float));
        }

        if (m_x_array.size() == 0) {
            cerr << endl
                 << "***Error! No particles found in binary file" << endl
                 << endl;
            throw runtime_error("Error extracting data from _binary file");
        }

        return;
    } else if (file_version == 3202) {
        f.read((char*)&m_input_position, sizeof(bool));
        f.read((char*)&m_input_velocity, sizeof(bool));
        //ct: add at 2014.03.08
        ReadExistTerm();

        int timestep;
        f.read((char*)&timestep, sizeof(unsigned int));
        m_timestep = timestep;
        unsigned int dimensions;
        f.read((char*)&dimensions, sizeof(unsigned int));
        m_num_dimensions = dimensions;
        float Lx, Ly, Lz;
        f.read((char*)&Lx, sizeof(float));
        f.read((char*)&Ly, sizeof(float));
        f.read((char*)&Lz, sizeof(float));
        m_box = vec(Lx, Ly, Lz);
        unsigned int np = 0;
        f.read((char*)&np, sizeof(unsigned int));
        if (m_only_parse_head)
            return;

        if (m_input_position) {
            m_x_array.resize(np);
            m_y_array.resize(np);
            m_z_array.resize(np);
            f.read((char*)&(m_x_array[0]), np * sizeof(float));
            f.read((char*)&(m_y_array[0]), np * sizeof(float));
            f.read((char*)&(m_z_array[0]), np * sizeof(float));
        }
        if (m_input_velocity) {
            m_vx_array.resize(np);
            m_vy_array.resize(np);
            m_vz_array.resize(np);
            f.read((char*)&(m_vx_array[0]), np * sizeof(float));
            f.read((char*)&(m_vy_array[0]), np * sizeof(float));
            f.read((char*)&(m_vz_array[0]), np * sizeof(float));
        }

        if (m_x_array.size() == 0) {
            cerr << endl
                 << "***Error! No particles found in binary file" << endl
                 << endl;
            throw runtime_error("Error extracting data from _binary file");
        }

        return;
    } else if (file_version == 1) {
        unsigned np = 0;
        f.read((char*)&m_input_position, sizeof(bool));
        f.read((char*)&m_input_velocity, sizeof(bool));
        f.read((char*)&m_input_image, sizeof(bool));
        f.read((char*)&m_input_type, sizeof(bool));
        f.read((char*)&m_input_charge, sizeof(bool));
        f.read((char*)&m_input_moleculeid, sizeof(bool));
        f.read((char*)&m_input_virial, sizeof(bool));
        f.read((char*)&m_input_force, sizeof(bool));
        f.read((char*)&m_input_anisotropy, sizeof(bool));
        f.read((char*)&m_input_init, sizeof(bool));
        f.read((char*)&m_input_cris, sizeof(bool));
        f.read((char*)&m_input_bond, sizeof(bool));
        f.read((char*)&m_input_angle, sizeof(bool));
        f.read((char*)&m_input_dihedral, sizeof(bool));
        f.read((char*)&m_input_integrator, sizeof(bool));

        ReadExistTerm();
        int timestep;
        f.read((char*)&timestep, sizeof(unsigned int));
        m_timestep = timestep;
        unsigned int dimensions;
        f.read((char*)&dimensions, sizeof(unsigned int));
        m_num_dimensions = dimensions;
        float Lx, Ly, Lz;
        f.read((char*)&Lx, sizeof(float));
        f.read((char*)&Ly, sizeof(float));
        f.read((char*)&Lz, sizeof(float));
        m_box = vec(Lx, Ly, Lz);
        f.read((char*)&np, sizeof(unsigned int));
        if (m_only_parse_head)
            return;

        m_tag_array.resize(np); // must be resized
        m_rtag_array.resize(np);
        f.read((char*)&(m_tag_array[0]), np * sizeof(unsigned int));
        f.read((char*)&(m_rtag_array[0]), np * sizeof(unsigned int));

        vector<vec4Dfloat> h_pos;
        vector<vec_int> h_image;
        vector<vec4Dfloat> h_vel;
        h_pos.resize(np);
        h_image.resize(np);
        h_vel.resize(np);

        if (m_input_position) {
            f.read((char*)&(h_pos[0]), 4 * np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_pos_array.push_back(vec(h_pos[tmpi].x, h_pos[tmpi].y, h_pos[tmpi].z));
            }
        }
        if (m_input_image) {
            f.read((char*)&(h_image[0]), 3 * np * sizeof(int));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_image_array.push_back(h_image[tmpi]);
            }
        }
        if (m_input_velocity) {
            f.read((char*)&(h_vel[0]), 4 * np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_vel_array.push_back(vec(h_vel[tmpi].x, h_vel[tmpi].y, h_vel[tmpi].z));
            }
        }

        return;
    } else if (file_version == 2) {

        unsigned np = 0;

        f.read((char*)&m_input_position, sizeof(bool));
        f.read((char*)&m_input_velocity, sizeof(bool));
        f.read((char*)&m_input_image, sizeof(bool));
        f.read((char*)&m_input_type, sizeof(bool));
        f.read((char*)&m_input_charge, sizeof(bool));
        f.read((char*)&m_input_body, sizeof(bool));
        f.read((char*)&m_input_virial, sizeof(bool));
        f.read((char*)&m_input_force, sizeof(bool));
        f.read((char*)&m_input_anisotropy, sizeof(bool));
        f.read((char*)&m_input_init, sizeof(bool));
        f.read((char*)&m_input_cris, sizeof(bool));
        f.read((char*)&m_input_bond, sizeof(bool));
        f.read((char*)&m_input_angle, sizeof(bool));
        f.read((char*)&m_input_dihedral, sizeof(bool));
        f.read((char*)&m_input_integrator, sizeof(bool));
        ReadExistTerm();

        int timestep;
        f.read((char*)&timestep, sizeof(unsigned int));
        m_timestep = timestep;
        unsigned int dimensions;
        f.read((char*)&dimensions, sizeof(unsigned int));
        m_num_dimensions = dimensions;
        float Lx, Ly, Lz;
        f.read((char*)&Lx, sizeof(float));
        f.read((char*)&Ly, sizeof(float));
        f.read((char*)&Lz, sizeof(float));
        m_box = vec(Lx, Ly, Lz);
        f.read((char*)&np, sizeof(unsigned int));
        if (m_only_parse_head)
            return;

        m_tag_array.resize(np); // must be resized
        m_rtag_array.resize(np);
        f.read((char*)&(m_tag_array[0]), np * sizeof(unsigned int));
        f.read((char*)&(m_rtag_array[0]), np * sizeof(unsigned int));

        vector<vec4Dfloat> h_pos;
        vector<vec_int> h_image;
        vector<vec4Dfloat> h_vel;
        h_pos.resize(np);
        h_image.resize(np);
        h_vel.resize(np);

        if (m_input_position) {
            f.read((char*)&(h_pos[0]), 4 * np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_pos_array.push_back(vec(h_pos[tmpi].x, h_pos[tmpi].y, h_pos[tmpi].z));
            }
        }
        if (m_input_image) {
            f.read((char*)&(h_image[0]), 3 * np * sizeof(int));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_image_array.push_back(h_image[tmpi]);
            }
        }
        if (m_input_velocity) {
            f.read((char*)&(h_vel[0]), 4 * np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_vel_array.push_back(vec(h_vel[tmpi].x, h_vel[tmpi].y, h_vel[tmpi].z));
            }
        }

        return;
    } else if (file_version == 3) {
        f.read((char*)&m_input_position, sizeof(bool));
        //cout<<"m_input_position"<<m_input_position<<endl;
        f.read((char*)&m_input_velocity, sizeof(bool));
        //cout<<"m_input_velocity"<<m_input_velocity<<endl;
        f.read((char*)&m_input_image, sizeof(bool));
        //cout<<"m_input_image"<<m_input_image<<endl;
        f.read((char*)&m_input_type, sizeof(bool));
        //cout<<"m_input_type"<<m_input_type<<endl;
        f.read((char*)&m_input_mass, sizeof(bool));
        //cout<<"m_input_mass"<<m_input_mass<<endl;
        f.read((char*)&m_input_diameter, sizeof(bool));
        //cout<<"m_input_diameter"<<m_input_diameter<<endl;
        f.read((char*)&m_input_body, sizeof(bool));
        //cout<<"m_input_body"<<m_input_body<<endl;
        f.read((char*)&m_input_accel, sizeof(bool));
        //cout<<"m_input_accel"<<m_input_accel<<endl;
        f.read((char*)&m_input_charge, sizeof(bool));
        //cout<<"m_input_charge"<<m_input_charge<<endl;
        f.read((char*)&m_input_moleculeid, sizeof(bool));
        //cout<<"m_input_moleculeid"<<m_input_moleculeid<<endl;
        f.read((char*)&m_input_virial, sizeof(bool));
        //cout<<"m_input_virial"<<m_input_virial<<endl;
        f.read((char*)&m_input_force, sizeof(bool));
        //cout<<"m_input_force"<<m_input_force<<endl;
        f.read((char*)&m_input_anisotropy, sizeof(bool));
        //cout<<"m_input_anisotropy"<<m_input_anisotropy<<endl;
        f.read((char*)&m_input_bond, sizeof(bool));
        //cout<<"m_input_bond"<<m_input_bond<<endl;
        f.read((char*)&m_input_angle, sizeof(bool));
        //cout<<"m_input_angle"<<m_input_angle<<endl;
        f.read((char*)&m_input_dihedral, sizeof(bool));
        //cout<<"m_input_dihedral"<<m_input_dihedral<<endl;
        f.read((char*)&m_input_integrator, sizeof(bool));
        //cout<<"m_input_integrator"<<m_input_integrator<<endl;
        f.read((char*)&m_input_rigid, sizeof(bool));
        //cout<<"m_input_rigid"<<m_input_rigid<<endl;

        //parse timestep
        ReadExistTerm();
        int timestep;
        f.read((char*)&timestep, sizeof(unsigned int));
        m_timestep = timestep;
        //cout<<"m_timestep"<<m_timestep<<endl;
        //parse dimensions
        unsigned int dimensions;
        f.read((char*)&dimensions, sizeof(unsigned int));
        m_num_dimensions = dimensions;
        //cout<<"m_num_dimensions"<<m_num_dimensions<<endl;
        //parse box
        float Lx, Ly, Lz;
        f.read((char*)&Lx, sizeof(float));
        f.read((char*)&Ly, sizeof(float));
        f.read((char*)&Lz, sizeof(float));
        m_box = vec(Lx, Ly, Lz);

        //allocate memory for particle arrays
        unsigned int np = 0;
        f.read((char*)&np, sizeof(unsigned int));
        if (m_only_parse_head)
            return;
        m_tag_array.resize(np);
        m_rtag_array.resize(np);
        //cout<<"np"<<np<<endl;
        //parse particle arrays
        f.read((char*)&(m_tag_array[0]), np * sizeof(unsigned int));
        f.read((char*)&(m_rtag_array[0]), np * sizeof(unsigned int));
        if (m_input_position) {
            m_x_array.resize(np);
            m_y_array.resize(np);
            m_z_array.resize(np);
            f.read((char*)&(m_x_array[0]), np * sizeof(float));
            f.read((char*)&(m_y_array[0]), np * sizeof(float));
            f.read((char*)&(m_z_array[0]), np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_pos_array.push_back(vec(m_x_array[tmpi], m_y_array[tmpi], m_z_array[tmpi]));
            }
        }
        if (m_input_image) {
            m_ix_array.resize(np);
            m_iy_array.resize(np);
            m_iz_array.resize(np);
            f.read((char*)&(m_ix_array[0]), np * sizeof(int));
            f.read((char*)&(m_iy_array[0]), np * sizeof(int));
            f.read((char*)&(m_iz_array[0]), np * sizeof(int));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_image_array.push_back(vec_int(m_ix_array[tmpi], m_iy_array[tmpi], m_iz_array[tmpi]));
            }
        }
        if (m_input_velocity) {
            m_vx_array.resize(np);
            m_vy_array.resize(np);
            m_vz_array.resize(np);
            f.read((char*)&(m_vx_array[0]), np * sizeof(float));
            f.read((char*)&(m_vy_array[0]), np * sizeof(float));
            f.read((char*)&(m_vz_array[0]), np * sizeof(float));
            for (unsigned i = 0; i != np; ++i) {
                unsigned tmpi = m_rtag_array[i];
                m_vel_array.push_back(vec(m_vx_array[tmpi], m_vy_array[tmpi], m_vz_array[tmpi]));
            }
        }
        if (m_input_accel) {
            m_ax_array.resize(np);
            m_ay_array.resize(np);
            m_az_array.resize(np);
            f.read((char*)&(m_ax_array[0]), np * sizeof(float));
            f.read((char*)&(m_ay_array[0]), np * sizeof(float));
            f.read((char*)&(m_az_array[0]), np * sizeof(float));
        }
        if (m_input_mass) {
            m_mass_array.resize(np);
            f.read((char*)&(m_mass_array[0]), np * sizeof(float));
        }
        if (m_input_diameter) {
            m_diameter_array.resize(np);
            f.read((char*)&(m_diameter_array[0]), np * sizeof(float));
        }
        if (m_input_charge) {
            m_charge_array.resize(np);
            f.read((char*)&(m_charge_array[0]), np * sizeof(float));
        }
        if (m_input_body) {
            m_body_array.resize(np);
            f.read((char*)&(m_body_array[0]), np * sizeof(unsigned int));
        }

        if (m_input_type) {
            m_idx_type_array.resize(np);
            //parse types
            unsigned int ntypes = 0;
            f.read((char*)&ntypes, sizeof(unsigned int));
            m_type_mapping.resize(ntypes);
            for (unsigned int i = 0; i < ntypes; i++)
                m_type_mapping[i] = read_string(f);
            f.read((char*)&(m_idx_type_array[0]), np * sizeof(unsigned int));
        }

        if (m_input_moleculeid) {
            m_Molecular_id.resize(np);
            //parse types
            unsigned int ntypes = 0;
            f.read((char*)&ntypes, sizeof(unsigned int));
            m_Molecular_mapping.resize(ntypes);
            for (unsigned int i = 0; i < ntypes; i++)
                m_Molecular_mapping[i] = read_string(f);
            f.read((char*)&(m_Molecular_id[0]), np * sizeof(unsigned int));
        }
        if (m_input_virial) {
            m_virial_array.resize(np);
            f.read((char*)&(m_virial_array[0]), np * sizeof(float));
        }

        if (m_input_force) {
            m_force_array.resize(np);
            f.read((char*)&(m_force_array[0]), 4 * np * sizeof(float));
        }
        if (m_input_anisotropy) {
            m_net_torque_array.resize(np);
            m_rotation_array.resize(np);
            m_orientation_array.resize(np);

            f.read((char*)&(m_net_torque_array[0]), 4 * np * sizeof(float));
            f.read((char*)&(m_rotation_array[0]), 4 * np * sizeof(float));
            f.read((char*)&(m_orientation_array[0]), 4 * np * sizeof(float));
        }

        //parse bonds
        if (m_input_bond) {
            unsigned int ntypes = 0;
            f.read((char*)&ntypes, sizeof(unsigned int));
            m_bond_type_mapping.resize(ntypes);
            for (unsigned int i = 0; i < ntypes; i++)
                m_bond_type_mapping[i] = read_string(f);

            /*ct:
          vector<unsigned int> tmp_bond_type_mapping;
          for(unsigned int i=0;i<m_bond_type_mapping.size();i++)
          {
          string::iterator iter=find(m_bond_type_mapping[i].begin(),m_bond_type_mapping[i].end(),'_');
          if(iter!=m_bond_type_mapping[i].end()) m_bond_type_mapping[i].erase(iter);
        //cout<<m_bond_type_mapping[i]<<endl;
        unsigned int a=lexical_cast<unsigned int>(m_bond_type_mapping[i]);
        tmp_bond_type_mapping.push_back(a);
        }
        //ct */

            unsigned int nb = 0;
            f.read((char*)&nb, sizeof(unsigned int));
            for (unsigned int j = 0; j < nb; j++) {
                unsigned int typ, a, b;
                f.read((char*)&typ, sizeof(unsigned int));
                f.read((char*)&a, sizeof(unsigned int));
                f.read((char*)&b, sizeof(unsigned int));

                m_bonds.push_back(Bond(m_bond_type_mapping[typ], a, b));
            }
        }
        //for(unsigned int i =0 ;i<m_bonds.size(); i++)
        //cout<<"m_bonds"<<m_bonds[i].a<<" "<<m_bonds[i].b<<endl;
        //parse angles
        if (m_input_angle) {
            unsigned int ntypes = 0;
            f.read((char*)&ntypes, sizeof(unsigned int));
            m_angle_type_mapping.resize(ntypes);
            for (unsigned int i = 0; i < ntypes; i++)
                m_angle_type_mapping[i] = read_string(f);
            /*ct:
          vector<unsigned int> tmp_angle_type_mapping;
          for(unsigned int i=0;i<m_angle_type_mapping.size();i++)
          {
          unsigned int a=lexical_cast<unsigned int>(m_angle_type_mapping[i]);
          tmp_angle_type_mapping.push_back(a);
          }
        //ct */

            unsigned int na = 0;
            f.read((char*)&na, sizeof(unsigned int));
            for (unsigned int j = 0; j < na; j++) {
                unsigned int typ, a, b, c;
                f.read((char*)&typ, sizeof(unsigned int));
                f.read((char*)&a, sizeof(unsigned int));
                f.read((char*)&b, sizeof(unsigned int));
                f.read((char*)&c, sizeof(unsigned int));
                m_angles.push_back(Angle(m_angle_type_mapping[typ], a, b, c));
            }
        }

        //parse dihedrals
        if (m_input_dihedral) {
            unsigned int ntypes = 0;
            f.read((char*)&ntypes, sizeof(unsigned int));
            m_dihedral_type_mapping.resize(ntypes);
            for (unsigned int i = 0; i < ntypes; i++)
                m_dihedral_type_mapping[i] = read_string(f);

            unsigned int nd = 0;
            f.read((char*)&nd, sizeof(unsigned int));
            for (unsigned int j = 0; j < nd; j++) {
                unsigned int typ, a, b, c, d;
                f.read((char*)&typ, sizeof(unsigned int));
                f.read((char*)&a, sizeof(unsigned int));
                f.read((char*)&b, sizeof(unsigned int));
                f.read((char*)&c, sizeof(unsigned int));
                f.read((char*)&d, sizeof(unsigned int));

                m_dihedrals.push_back(Dihedral(typ, a, b, c, d));
            }
        }

        //parse integrator states
        if (m_input_integrator) {
            std::vector<IntegratorVariables> v;
            unsigned int ni = 0;
            f.read((char*)&ni, sizeof(unsigned int));
            v.resize(ni);
            for (unsigned int j = 0; j < ni; j++) {
                v[j].type = read_string(f);

                v[j].variable.clear();
                unsigned int nv = 0;
                f.read((char*)&nv, sizeof(unsigned int));
                for (unsigned int k = 0; k < nv; k++) {
                    float var;
                    f.read((char*)&var, sizeof(float));
                    v[j].variable.push_back(var);
                }
            }
            m_integrator_variables = v;
        }

        // parse rigid bodies
        if (m_input_rigid) {
            unsigned int n_bodies = 0;
            f.read((char*)&n_bodies, sizeof(unsigned int));

            if (n_bodies == 0)
                return;

            m_com.resize(n_bodies);
            m_vel.resize(n_bodies);
            m_angmom.resize(n_bodies);
            m_body_image.resize(n_bodies);

            for (unsigned int body = 0; body < n_bodies; body++) {
                f.read((char*)&(m_com[body].x), sizeof(float));
                f.read((char*)&(m_com[body].y), sizeof(float));
                f.read((char*)&(m_com[body].z), sizeof(float));
                f.read((char*)&(m_com[body].w), sizeof(float));

                f.read((char*)&(m_vel[body].x), sizeof(float));
                f.read((char*)&(m_vel[body].y), sizeof(float));
                f.read((char*)&(m_vel[body].z), sizeof(float));
                f.read((char*)&(m_vel[body].w), sizeof(float));

                f.read((char*)&(m_angmom[body].x), sizeof(float));
                f.read((char*)&(m_angmom[body].y), sizeof(float));
                f.read((char*)&(m_angmom[body].z), sizeof(float));
                f.read((char*)&(m_angmom[body].w), sizeof(float));

                f.read((char*)&(m_body_image[body].x), sizeof(int));
                f.read((char*)&(m_body_image[body].y), sizeof(int));
                f.read((char*)&(m_body_image[body].z), sizeof(int));
            }
        }

        // check for required items in the file
        if (m_x_array.size() == 0) {
            cerr << endl
                 << "***Error! No particles found in binary file" << endl
                 << endl;
            throw runtime_error("Error extracting data from _binary file");
        }
    } else {
        throw runtime_error("Unknown version of binary file.");
    }
}

unsigned int BinReader::getNumBondTypes() const
{
    return (unsigned int)m_bond_type_mapping.size();
}
unsigned int BinReader::getNumAngleTypes() const
{
    return (unsigned int)m_angle_type_mapping.size();
}
unsigned int BinReader::getNumDihedralTypes() const
{
    return (unsigned int)m_dihedral_type_mapping.size();
}

int BinReader::getFileVersion() const
{
    return file_version;
}

void export_BinReader()
{
    class_<BinReader, bases<Reader>, boost::noncopyable>("BinReader", init<const string&>())
        .def(init<const string&, bool>())
        .def("getFileVersion", &BinReader::getFileVersion);
}
