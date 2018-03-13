/**************************************************************************
 *   CMakeLists.txt  --  This file is part of lanius.                     *
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   lanius is free software: you can redistribute it and/or modify       *
 *   it under the terms of the GNU General Public License as published    *
 *   by the Free Software Foundation, either version 3 of the License,    *
 *   or (at your option) any later version.                               *
 *                                                                        *
 *   lanius is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include <iostream>
#include <chrono>
#include <tclap/CmdLine.h>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"
#include "surface_creator.h"
#include "crystal_database.h"

int main(int argc, char *argv[]) {

    try {
        TCLAP::CmdLine cmd("Create surface given unit cell and miller indices", ' ', PROGRAM_VERSION);

        //**************************************
        // Inform user about execution
        //**************************************
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing Lanius v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        // keep track of execution time
        auto start = std::chrono::system_clock::now();

        // output filename
        TCLAP::ValueArg<std::string> arg_output_filename("o","filename","Filename to coordinates to",true,"POSCAR","string");
        cmd.add(arg_output_filename);

        // miller_indices
        TCLAP::ValueArg<std::string> arg_m("m","miller","Miller indices",true, "1,1,1","string");
        cmd.add(arg_m);

        // crystal type
        TCLAP::ValueArg<std::string> arg_c("c","crystal","Crystal structure",true, "fcc","string");
        cmd.add(arg_c);

        // element
        TCLAP::ValueArg<std::string> arg_e("e","element","Element type",true, "Co","string");
        cmd.add(arg_e);

        // dimensions
        TCLAP::ValueArg<std::string> arg_d("d","dimensions","Dimension of surface",true, "2,2,4","string");
        cmd.add(arg_d);

        // vacuum
        TCLAP::ValueArg<std::string> arg_v("v","vacuum","Vacuum thickness in A", true, "10","string");
        cmd.add(arg_v);

        cmd.parse(argc, argv);

        //**************************************
        // parsing values
        //**************************************
        const std::string filename = arg_output_filename.getValue();
        const std::string crystal = arg_c.getValue();
        const double vacuum = boost::lexical_cast<double>(arg_v.getValue());
        const std::string element = arg_e.getValue();

        //**************************************
        // parsing indices
        //**************************************
        const boost::regex re_scalar_triplet("^([0-9]+),([0-9]+),([0-9]+)$");
        std::vector<int> dimensions(6,0);
        boost::smatch what;
        const std::string miller = arg_m.getValue();
        if(boost::regex_match(miller, what, re_scalar_triplet)) {
            dimensions[0] = boost::lexical_cast<int>(what[1]);
            dimensions[1] = boost::lexical_cast<int>(what[2]);
            dimensions[2] = boost::lexical_cast<int>(what[3]);
        }
        const std::string dims = arg_d.getValue();
        if(boost::regex_match(dims, what, re_scalar_triplet)) {
            dimensions[3] = boost::lexical_cast<int>(what[1]);
            dimensions[4] = boost::lexical_cast<int>(what[2]);
            dimensions[5] = boost::lexical_cast<int>(what[3]);
        }

        std::cout << "Constructing " << crystal << " lattice." << std::endl;
        std::cout << "Cutting with surface normal: (" << dimensions[0] << "," << dimensions[1] << "," << dimensions[2] << ")." << std::endl;
        std::cout << "Expanding surface with: (" << dimensions[3] << "," << dimensions[4] << "," << dimensions[5] << ")." << std::endl;
        std::cout << "Building vacuum layer: " << vacuum << " angstrom." << std::endl;

        SurfaceCreator sf;
        CrystalDatabase cdb;
        const std::vector<double> lc = cdb.get_lattice_parameters(element + "-" + crystal);
        sf.construct_unitcell(crystal, element, lc[0], lc[1], lc[2]);
        sf.cleave(dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4], dimensions[5]);
        sf.create_surface(vacuum);
        sf.export_surface(filename);

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Done in " << elapsed_seconds.count() << " seconds" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
