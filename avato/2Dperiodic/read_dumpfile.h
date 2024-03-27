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
# ifndef READ_DUMPFILE_H
# define READ_DUMPFILE_H

# include <cassert>
# include <cstdio>
# include <cstring>
# include <cstdlib>

# ifndef MAXLINELEN
#   define MAXLINELEN 65535
# endif

# ifndef NDIM
#   define NDIM 2
#   define NDIM_UNDEF
# endif
# ifndef RANDOMMOTORS
#   define RANDOMMOTORS 0
#   define RANDOMMOTORS_UNDEF
# endif
# ifndef VARISPHERES
#   define VARISPHERES 0
#   define VARISPHERES_UNDEF
# endif
# ifndef BPARTICLE
#   define BPARTICLE 0
#   define BPARTICLE_UNDEF
# endif

template<class Box, class Point>
bool read_dumpfile(Box& domain, std::vector<Point>& coms, std::vector<double>& radii, const char* dumpfilename) // <<<
{
  // open file or pipe <<<
  FILE* fp = NULL;
  bool usepipe = (strncmp(dumpfilename + strlen(dumpfilename) - 3, ".gz", 3) == 0);
  if (usepipe) {
    char command[MAXLINELEN];
    // What happens if NFS hangs? Test first whether file is readable.
    fp = fopen(dumpfilename, "rb");
    if (fp == NULL)
        return false;
    fclose(fp);
    // then launch gunzip
    sprintf(command, "gunzip -c '%s'", dumpfilename);
    fp = popen(command, "r");
  } else {
    fp = fopen(dumpfilename, "r");
  }
  if (fp == NULL) {
    //fprintf(stderr, "Cannot open file \"%s\".\n", dumpfilename);
    return false;
  } // >>>

  // read header <<<
  char line_array[MAXLINELEN];
  char* line = line_array;
  int nStdSpheres=0;
  //int nInitGrid_x, nInitGrid_y;
  double Period1_x, Period1_y, Period2_x, Period2_y;
  # if NDIM == 3
  //int nInitGrid_z;
  double Period1_z, Period2_z, Period3_x, Period3_y, Period3_z;
  # endif
  bool haveveloc = false;
  bool haveperiod = false;
  while (strncmp(fgets(line, MAXLINELEN, fp), "ITEM: ATOMS", 11) != 0) {
    if (strncmp(line, "ITEM: TIMESTEP", 14) == 0) {
      fgets(line, MAXLINELEN, fp);
    } else if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) == 0) {
      if (fscanf(fp, "%d\n", &nStdSpheres) != 1) {fprintf(stderr, "ERROR at NUMBER OF ATOMS\n"); exit(1);}
    } else if (strncmp(line, "ITEM: TEMPERATURE", 17) == 0) {
      fgets(line, MAXLINELEN, fp);
      //if (fscanf(fp, "%lg\n", &kT) != 1) {fprintf(stderr, "ERROR at TEMPERATURE\n"); exit(1);}
    } else if (strncmp(line, "ITEM: INITIAL GRID\n", MAXLINELEN) == 0) {
      fgets(line, MAXLINELEN, fp);
    } else if (strncmp(line, "ITEM: PERIOD VECTORS\n", MAXLINELEN) == 0) {
      // XXX we use this field for the *current* grid vectors in the case of ONLINEDEFORM
      # if NDIM == 3
      if (fscanf(fp, "%lg %lg %lg\n%lg %lg %lg\n%lg %lg %lg\n",
                 &(Period1_x), &(Period1_y), &(Period1_z),
                 &(Period2_x), &(Period2_y), &(Period2_z),
                 &(Period3_x), &(Period3_y), &(Period3_z)) != 9) {fprintf(stderr, "ERROR at PERIOD VECTORS\n"); exit(1);}
      // period vectors must be orthogonal in CGAL
      assert(Period1_y == 0.0 and Period1_z == 0.0);
      assert(Period2_x == 0.0 and Period2_z == 0.0);
      assert(Period3_x == 0.0 and Period3_y == 0.0);
      # else // NDIM == 2
      if (fscanf(fp, "%lg %lg\n%lg %lg\n",
                 &(Period1_x), &(Period1_y),
                 &(Period2_x), &(Period2_y)) != 4) {fprintf(stderr, "ERROR at PERIOD VECTORS\n"); exit(1);}
      assert(Period1_y == 0.0);
      assert(Period2_x == 0.0);
      # endif
    } else if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      for (int i=0; i<NDIM; i++) fgets(line, MAXLINELEN, fp);
    } else if (strncmp(line, "ITEM: CONTENT", 13) == 0) {
      fgets(line, MAXLINELEN, fp);
      assert(strncmp(line, "id type com", 11) == 0);
      line += 11;
      # if RANDOMMOTORS
      assert(strncmp(line, " ori", 4) == 0);
      line += 4;
      # endif
      haveveloc = (strncmp(line, " vel", 4) == 0);
      if (haveveloc) line += 4;
      haveperiod = (strncmp(line, " periods", 8) == 0);
      if (haveperiod) line += 8;
      # if VARISPHERES or BPARTICLE
      assert(strncmp(line, " mass radius", 12) == 0);
      # endif
    } else if (strncmp(line, "ITEM: QUALIFIERS", 16) == 0) {
      fgets(line, MAXLINELEN, fp);
    } else {
      fprintf(stderr, "ERROR: unknown line:\n%s\n", line);
      exit(EXIT_FAILURE);
    }
  } // >>>

  // create a new periodic domain: must be orthogonal in CGAL
  # if NDIM == 3
  domain = Box(-0.5*Period1_x, -0.5*Period2_y, -0.5*Period3_z,
                0.5*Period1_x,  0.5*Period2_y,  0.5*Period3_z);
  # elif NDIM == 2
  domain = Box(-0.5*Period1_x, -0.5*Period2_y,
                0.5*Period1_x,  0.5*Period2_y);
  # endif

  // read the positions + radii <<<
  assert(not haveveloc);
  coms.resize(nStdSpheres);
  radii.resize(nStdSpheres);
  for (int i=0; i<nStdSpheres; i++) {
    fgets(line, MAXLINELEN, fp);
    int cnt=-1, type=-1;
    double mass, radius;
    double com_x, com_y;

    #define ERRR(str) \
    { \
      fprintf(stderr, "%s\n%sn", str, line); \
      exit(EXIT_FAILURE); \
    }

    # if NDIM == 3
    double com_z;
    if (haveperiod) {
      int px, py, pz;
      if (sscanf(line, "%d %d %lg %lg %lg %d %d %d %lg %lg\n", &cnt, &type, &com_x, &com_y, &com_z, &px, &py, &pz, &mass, &radius) != 10) ERRR("Reading failed");
    } else {
      if (sscanf(line, "%d %d %lg %lg %lg %lg %lg\n", &cnt, &type, &com_x, &com_y, &com_z, &mass, &radius) != 7) ERRR("Reading failed");
    }
    # else
    if (haveperiod) {
      int px, py;
      if (sscanf(line, "%d %d %lg %lg %d %d %lg %lg\n", &cnt, &type, &com_x, &com_y, &px, &py, &mass, &radius) != 8) ERRR("Reading failed");
    } else {
      if (sscanf(line, "%d %d %lg %lg %lg %lg\n", &cnt, &type, &com_x, &com_y, &mass, &radius) != 6) ERRR("Reading failed");
    }
    # endif

    # if NDIM == 3
    coms[i] = Point(com_x, com_y, com_z);
    # elif NDIM == 2
    coms[i] = Point(com_x, com_y);
    # endif
    radii[i] = radius;
  } // >>>

  if (usepipe) pclose(fp);
  else fclose(fp);

  return true;
} // >>>

# ifdef NDIM_UNDEF
#   undef NDIM
# endif
# ifdef RANDOMMOTORS_UNDEF
#   undef RANDOMMOTORS
# endif
# ifdef VARISPHERES_UNDEF
#   undef VARISPHERES
# endif
# ifdef BPARTICLE_UNDEF
#   undef BPARTICLE
# endif
# endif // READ_DUMPFILE_H

// vim:foldmethod=marker:foldmarker=<<<,>>>
