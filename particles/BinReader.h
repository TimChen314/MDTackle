#include "MolDataStruct.h"
#include "Reader.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;

#ifndef __BIN_READER_H__
#define __BIN_READER_H__

/* ct: The points may need to promote:
   1. code a Template function to
*/

class BinReader : public Reader {
public:
    BinReader(const std::string& fname);
    BinReader(const std::string& fname, bool m_only_parse_head);
    BinReader(const std::string& fname, const unsigned block_begin_, const unsigned block_end_);

    virtual unsigned int getNumParticles() const;
    virtual unsigned int getNumParticleTypes() const;
    virtual vec getBox() const;
    virtual vector<vec>& getRealPos();
    virtual vector<vec>& getBoxPos();

    //! Returns the timestep of the simulation
    virtual unsigned int getTimeStep() const;

    //! Sets the timestep of the simulation
    virtual void setTimeStep(unsigned int ts);

    //! Returns the number of dimensions
    virtual unsigned int getNumDimensions() const;

    //! Returns the number of bond types to be created
    virtual unsigned int getNumBondTypes() const;

    //! Returns the number of angle types to be created
    virtual unsigned int getNumAngleTypes() const;

    //! Returns the number of dihedral types to be created
    virtual unsigned int getNumDihedralTypes() const;

    //! Returns the file_version
    virtual int getFileVersion() const;

private:
    //! Helper function to read the input file
    void readFile(const std::string& fname);
    void convertToVec();
    void ReadExistTerm();

    //ct: Special array for BinReader, note that these arrays may not exist:
    //ct: since m_ax_array's written as float type, it can only be read by float type. This has been validated by myself.
    std::vector<unsigned int> m_tag_array;  //!< tags of all particles loaded
    std::vector<unsigned int> m_rtag_array; //!< inverse tags of all particles loaded
    std::vector<float> m_ax_array;          //!< x acceleration of all particles loaded
    std::vector<float> m_ay_array;          //!< y acceleration of all particles loaded
    std::vector<float> m_az_array;          //!< z acceleration of all particles loaded
    std::vector<float> m_x_array;           //!< x position of all particles loaded
    std::vector<float> m_y_array;           //!< y position of all particles loaded
    std::vector<float> m_z_array;           //!< z position of all particles loaded
    std::vector<float> m_vx_array;          //!< x velocity of all particles loaded
    std::vector<float> m_vy_array;          //!< y velocity of all particles loaded
    std::vector<float> m_vz_array;          //!< z velocity of all particles loaded
    std::vector<int> m_ix_array;            //!< x image of all particles loaded
    std::vector<int> m_iy_array;            //!< y image of all particles loaded
    std::vector<int> m_iz_array;            //!< z image of all particles loaded

    int file_version;                                        //ct:
    std::vector<IntegratorVariables> m_integrator_variables; //!< Integrator variables read in from file

    std::vector<vec4Dfloat> m_com;     //!< n_bodies length 1D array of center of mass positions
    std::vector<vec4Dfloat> m_vel;     //!< n_bodies length 1D array of body velocities
    std::vector<vec4Dfloat> m_angmom;  //!< n_bodies length 1D array of angular momenta in the space frame
    std::vector<vec_int> m_body_image; //!< n_bodies length 1D array of the body image

    bool m_input_position = false; //!< true if the particle positions should be written
    bool m_input_image = false;
    bool m_input_velocity = false;
    bool m_input_type = false;
    bool m_input_mass = false;
    bool m_input_diameter = false;
    bool m_input_body = false;
    bool m_input_accel = false;
    bool m_input_charge = false;
    bool m_input_moleculeid = false;
    bool m_input_virial = false;
    bool m_input_force = false;
    bool m_input_anisotropy = false;
    bool m_input_bond = false;
    bool m_input_angle = false;
    bool m_input_dihedral = false;
    bool m_input_integrator = false;
    bool m_input_rigid = false;
    bool m_input_cris = false;
    bool m_input_init = false;
    bool m_only_parse_head = false; //ct:  if you don't need read position and other arrays, then set this bool varible true.
    bool m_read_by_block = false;

    unsigned block_begin, block_end; //ct: used for reading by block
};

void export_BinReader();

#endif
