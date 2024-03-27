// Copyright (c) 2015 CNRS (France).
// Author: Michael Schindler <michael.schindler@espci.fr>
// All rights reserved.
//
// This file is part of AVATO (https://www.pct.espci.fr/~michael/avato)
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <csignal>
#include <sys/stat.h>
#include "Periodic_2_dual_triangulations_2.h"
#include <libgen.h>
#include <vector>
#include <random>

//pedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef typename K::FT FT;
typedef Periodic_2_dual_triangulations_2<K> TriVor;

bool terminate = false;
void TERMhandler(__attribute__((unused)) int sig) // <<<
// handler for SIGTERM, -15
{
  terminate = true;
} // >>>

int main(int argc, char* argv[])
{
  FT classes_prec = -1.0; // default: do not make radii classes
  std::string myname(basename(argv[0]));
  // artificial treatment of radii? <<<
  std::vector<double> overwrite_radii;
  bool mono = (myname.find_last_of("_") != static_cast<size_t>(-1) and myname.substr(myname.find_last_of("_"), myname.back()) == "_mono");
  if (mono)
  {
    std::cout << "Pretending that all radii are equal (monodisperse)." << std::endl;
    classes_prec = 1.0e-3; // there will be only one class. Determine it automatically
  }

  bool shuffle = (myname.find_last_of("_") != static_cast<size_t>(-1) and myname.substr(myname.find_last_of("_"), myname.back()) == "_shuffle");
  if (shuffle)
  {
    std::cout << "Shuffling all radii." << std::endl;
  }
  // >>>

  TriVor TV;
  bool keep_domain = false;
  signal(SIGTERM, TERMhandler);
  signal(SIGUSR1, TERMhandler);
  for (int i=1; i<argc; i++)
  {
    std::string dumpname(argv[i]);
    std::string outname = dumpname.substr(0, dumpname.find_last_of(".")) + ".av" + (mono ? "_mono" : "") + (shuffle ? "_shuffle" : "");
    struct stat s;
    bool skip_this = false;
    if (stat(outname.c_str(), &s) == 0)
    {
      std::cout << "WARNING: File " << outname << " exists already. Skipping treatment of " << argv[i] << std::endl;
      skip_this = true;
    }
    else if (stat((outname + ".gz").c_str(), &s) == 0)
    {
      std::cout << "WARNING: File " << outname << ".gz exists already. Skipping treatment of " << argv[i] << std::endl;
      skip_this = true;
    }

    // XXX even if we skip some files, read radii and geometry from the first!
    if (skip_this and i > 1)
      continue;

    if (not TV.read_points(argv[i], keep_domain)) {
        std::cout << "Cannot read \"" << argv[i] << "\"." << std::endl;
        continue;
    }
    std::cout << "Treating \"" << argv[i] << "\" ..." << std::endl;
    std::cout << "  bb_min = " << TV.domain().xmin() << " " << TV.domain().ymin() << std::endl
              << "  bb_max = " << TV.domain().xmax() << " " << TV.domain().ymax() << std::endl
              << "  ndens = " << static_cast<int>(TV.orig_points().size()) / ((TV.domain().xmax()-TV.domain().xmin()) * (TV.domain().ymax()-TV.domain().ymin())) << std::endl;

    // init the overwrite_radii <<<
    if (i == 1) // first iteration
    {
      if (mono) {
        overwrite_radii = TV.orig_radii();
        //// set average radius everywhere
        //double radius = 0.0;
        //for (size_t i=0; i<overwrite_radii.size(); i++)
        //  radius += overwrite_radii[i];
        //radius /= overwrite_radii.size();
        double radius = 0.5;
        std::cout << "INFO: Using radius = " << radius << " for pseudo-monodisperse." << std::endl;
        overwrite_radii.assign(overwrite_radii.size(), radius);
      }
      if (shuffle) {
        overwrite_radii = TV.orig_radii();
        //std::default_random_engine rand;
        std::mt19937_64 rand;
        rand.seed(43); // XXX use same seed because we often interrupt series of dump.gz -> dump.av_shuffle
        std::uniform_int_distribution<size_t> uniform(0, overwrite_radii.size()-1);
        for (int repeats=0; repeats<2; repeats++)
        {
          for (size_t i=0; i<overwrite_radii.size(); i++)
          {
            // swap items i and j:
            size_t j = uniform(rand);
            double tmp = overwrite_radii[i];
            overwrite_radii[i] = overwrite_radii[j];
            overwrite_radii[j] = tmp;
          }
        }
      }
    }
    // >>>

    if (skip_this)
      continue;

    if (mono or shuffle) {
      std::cout << "INFO: Overwriting the radii." << std::endl;
      assert(overwrite_radii.size() == TV.orig_radii().size());
      TV.overwrite_radii(overwrite_radii);
    }
    std::cout << "INFO: Writing data to " << outname << std::endl;
    std::ofstream outfile;
    outfile.open(outname);
    TV.treat_available_volumes(true, outfile, classes_prec);
    outfile.close();

    keep_domain = true;

    if (terminate) break;
  }
  return 0;
}

// vim:foldmethod=marker:foldmarker=<<<,>>>
