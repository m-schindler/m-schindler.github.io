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
#include "Periodic_2_dual_triangulations_2.h"

//pedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef typename K::FT FT;
typedef Periodic_2_dual_triangulations_2<K> TriVor;

int main()
{
  TriVor TV;
  TV.unit_tests();
  return 0;
}

// vim:foldmethod=marker:foldmarker=<<<,>>>
