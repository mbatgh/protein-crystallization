
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include "aspbe.hpp"

void read_parameters(const int argc, const char** argv, s_pars& my_pars)
{
    namespace po = boost::program_options;

    po::options_description generic("Command line options");
    generic.add_options()
      ("version,v", "print version")
      ("help,h", "produce help message")

      ("confil,f",        po::value< std::string >( &my_pars.confil ),  "configuration file name" )
      ("runid,i",         po::value< std::string >( &my_pars.runid ),   "run-id, for naming of output files" )
      ("rseed,r",         po::value< int    >(&my_pars.rseed ),         "seed for random number generator" )
      ("fr_AS,a",         po::value< double >(&my_pars.fr_AS ),         "flow rate of antisolvent [cm3/sec]" )
      ("na",              po::value< double >(&my_pars.k_nucl_a ),      "nucleation na,  B = na ((c-csat)/csat)^nb [#/cm3/sec]" )
      ("ga",              po::value< double >(&my_pars.k_growth_a ),    "growth ga,      G = ga ((c-csat)/csat)^gb [cm/sec]" )
      ("nb",              po::value< double >(&my_pars.k_nucl_b ),      "nucleation nb,  B = na ((c-csat)/csat)^nb [#/cm3/sec]" )
      ("gb",              po::value< double >(&my_pars.k_growth_b ),    "growth gb,      G = ga ((c-csat)/csat)^gb [cm/sec]" )
      ("ak",              po::value< double >(&my_pars.ak ),            "aggregation kernel const" )
      ;

    po::options_description config("Various options");
    config.add_options()
      ("name",            po::value< std::string >( &my_pars.name ),    "protein name" )
      ("asname",          po::value< std::string >( &my_pars.asname ),  "anti-solvent name" )
      ("rseed",           po::value< int    >(&my_pars.rseed ),         "seed for random number generator" )
      ("DOCLD",           po::value< int    >(&my_pars.DOCLD ),         "1: generate CLD, 0: nope" )
      ("sizethresh",      po::value< double >(&my_pars.sizethresh),     "the desired minimum crystal size [cm]" )
      ("nrot",            po::value< int    >(&my_pars.nrot ),          "number of random rotations for CSD2CLD" )
      ("SizeNucl",        po::value< double >(&my_pars.SizeNucl ),      "average size of nucleus [cm]" )
      ("ga",              po::value< double >(&my_pars.k_growth_a ),    "growth kg, G = kg ((c-csat)/csat)^dg [cm/sec]" )
      ("gb",              po::value< double >(&my_pars.k_growth_b ),    "growth dg, G = kg ((c-csat)/csat)^dg []" )
      ("na",              po::value< double >(&my_pars.k_nucl_a ),      "nucleation kn, G = kn ((c-csat)/csat)^dn [#/cm3/sec]" )
      ("nb",              po::value< double >(&my_pars.k_nucl_b ),      "nucleation dn, G = kn ((c-csat)/csat)^dn []" )
      ("pcp0",            po::value< double >(&my_pars.pcp0),           "percentage of inital protein conc that determines tpc0" )
      ("ds_protein",      po::value< double >(&my_pars.ds_protein ),    "density protein in crystal [g/cm3]" )
      ("ds_water",        po::value< double >(&my_pars.ds_water),       "density water in crystal [g/cm3]" )
      ("printlevel",      po::value< int >   (&my_pars.printlevel),     "boolean, 0: print nothing, 1: print CSD of tn, 2: print CSD of all t" )
      ("dt_0",            po::value< double >(&my_pars.dt_0 ),          "initial time step [sec]" )
      ("dt_write",        po::value< double >(&my_pars.dt_write ),      "write interval for running output [sec]" )
      ("t_end",           po::value< double >(&my_pars.t_end ),         "total simulation time [sec]" )
      ("sizeexp",         po::value< double >(&my_pars.sizeexp ),       "exponent for grid size generation" )
      ("typecsd0",        po::value< int >(   &my_pars.typecsd0 ),      "type of seed CSD (0: none, 1: square, 2: normal dist)" )
      ("mcsd0",           po::value< double >(&my_pars.mcsd0 ),         "average size of seed CSD" )
      ("scsd0",           po::value< double >(&my_pars.scsd0 ),         "stdev/halfwidth of seed CSD" )
      ("ms_protein_0",    po::value< double >(&my_pars.ms_protein_0 ),  "mass of protein in seed" )
      ("shapefac",        po::value< double >(&my_pars.shapefac ),      "shapefactor ratio (crystal volume)/cube" )
      ("fr_AS",           po::value< double >(&my_pars.fr_AS ),         "flow rate of antisolvent [cm3/sec]" )
      ("fr_AS2",          po::value< double >(&my_pars.fr_AS2 ),        "flow rate of antisolvent after t1 [cm3/sec]" )
      ("fr_AS3",          po::value< double >(&my_pars.fr_AS3 ),        "flow rate of antisolvent agter t2 [cm3/sec]" )
      ("t1",              po::value< double >(&my_pars.t1 ),            "time to change to fr_AS2 [sec]" )
      ("t2",              po::value< double >(&my_pars.t2 ),            "time to change to fr_AS3 [sec]" )
      ("Temp",            po::value< double >(&my_pars.Temp ),          "Temperature [degC]" )
      ("pH",              po::value< double >(&my_pars.pH ),            "pH []" )
      ("ml_protein_0",    po::value< double >(&my_pars.ml_protein_0 ),  "mass of protein in solution at t=0 [g]" )
      ("ml_AS_max",       po::value< double >(&my_pars.ml_AS_max ),     "total volume of antisolvent to add" )
      ("ml_water_0",      po::value< double >(&my_pars.ml_water_0 ),    "mass of water in solution at t=0 [g]" )
      ("nclasses",        po::value< int    >(&my_pars.NClasses_max ),  "number of grid points" )
      ("dl_water_res",    po::value< double >(&my_pars.dl_water_res ),  "density of water in AS solution [g/cm3]" )
      ("dl_AS_res",       po::value< double >(&my_pars.dl_AS_res ),     "density of water in AS solution [g/cm3]" )
      ("const_supsat",    po::value< double >(&my_pars.const_supsat ),  "constant super saturation (for testing)" )
      ("at",              po::value< int    >(&my_pars.at ),            "aggregation type - 0: none, 1: const, 2: sum, 3: product, 4: (Ri+Rj)^3, 5: G*(Ri+Rj)^3")
      ("ak",              po::value< double >(&my_pars.ak ),            "aggregation kernel const" )
      ;
  
    po::options_description cmdline_options;
    cmdline_options.add(generic);
    po::options_description various_options;
    various_options.add(config);

    po::variables_map vm;

    store( po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    notify(vm);

    if (vm.count("help")) {
      cout << generic << endl;
      cout << config << endl;
      exit(1);
    }

    if (vm.count("version")) {
      cout << "aspbe version 2.01" << endl;
      exit(1);
    }

    if (!vm.count("runid")) {
      cout << endl << "no runid given!" << endl << endl;
      cout << generic << endl;
      cout << config << endl;
      exit(0);
    }

    ifstream ifs(my_pars.confil.c_str());
    if (!ifs)
    {
        cout << "can not open config file: " << my_pars.confil << "\n";
        exit(0);
    }
    else
    {
        store(parse_config_file(ifs, various_options), vm);
        notify(vm);
    }
}
