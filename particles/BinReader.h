#include "MolDataStruct.h"
#include "Reader.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;

#ifndef __BIN_READER_H__
#define __BIN_READER_H__


/* ct: The points may need to promote:
   1. code a Template function to
*/

class BinReader : public Reader
    {
    public:

        BinReader(const std::string &fname);
        BinReader(const std::string &fname,bool m_only_parse_head);

        virtual unsigned int getNumParticles() const;
        virtual unsigned int getNumParticleTypes() const;
        virtual vec getBox() const;
  virtual vector< vec >& getRealPos() ;
  virtual vector<vec >& getBoxPos() ;



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
        void readFile(const std::string &fname);
        void convertToVec();
        void ReadExistTerm();

//ct: Special array for BinReader, note that these arrays may not exist:
//ct: since m_ax_array's written as float type, it can only be read by float type. This has been validated by myself.
        std::vector< unsigned int > m_tag_array;     //!< tags of all particles loaded
        std::vector< unsigned int > m_rtag_array;    //!< inverse tags of all particles loaded
        std::vector< float > m_ax_array;           //!< x acceleration of all particles loaded
        std::vector< float > m_ay_array;           //!< y acceleration of all particles loaded
        std::vector< float > m_az_array;           //!< z acceleration of all particles loaded
        std::vector< float > m_x_array;            //!< x position of all particles loaded
        std::vector< float > m_y_array;            //!< y position of all particles loaded
        std::vector< float > m_z_array;            //!< z position of all particles loaded
        std::vector< float > m_vx_array;           //!< x velocity of all particles loaded
        std::vector< float > m_vy_array;           //!< y velocity of all particles loaded
        std::vector< float > m_vz_array;           //!< z velocity of all particles loaded
        std::vector< int > m_ix_array;              //!< x image of all particles loaded
        std::vector< int > m_iy_array;              //!< y image of all particles loaded
        std::vector< int > m_iz_array;              //!< z image of all particles loaded


  int file_version; //ct:
  std::vector<IntegratorVariables> m_integrator_variables; //!< Integrator variables read in from file


        std::vector< vec4Dfloat> m_com;                    //!< n_bodies length 1D array of center of mass positions
        std::vector< vec4Dfloat> m_vel;                    //!< n_bodies length 1D array of body velocities
        std::vector< vec4Dfloat> m_angmom;                 //!< n_bodies length 1D array of angular momenta in the space frame
        std::vector< vec_int > m_body_image;                //!< n_bodies length 1D array of the body image

        bool m_imput_position;     //!< true if the particle positions should be written
        bool m_imput_image;        //!< true if the particle images should be written
        bool m_imput_velocity;     //!< true if the particle velocities should be written
        bool m_imput_type;         //!< true if the particle types should be written
        bool m_imput_mass;         //!< true if the particle masses should be written
        bool m_imput_diameter;     //!< true if the particle diameters should be written
	bool m_imput_body;			//!< true if the particle body should be written
        bool m_imput_accel;        //!< true if acceleration should be written
        bool m_imput_charge;       //!< true if charge should be written
        bool m_imput_moleculeid;	//!< true if particle's molecular id should be written
        bool m_imput_virial;       //!< true if virial should be written
        bool m_imput_force;        //!< true if force should be written
        bool m_imput_anisotropy;  //!< true if orientation should be written
        bool m_imput_bond;         //!< true if the bonds should be writte
        bool m_imput_angle;        //!< true if the angles should be written
        bool m_imput_dihedral;     //!< true if dihedrals should be written
        bool m_imput_integrator;	//!< true if integrators should be written
        bool m_imput_rigid;		//!< true if rigids should be written
        bool m_imput_cris;
        bool m_imput_init;
        bool m_only_parse_head;//ct:  if you don't need read position and other arrays, then set this bool varible true.

    };

void export_BinReader();


#endif



