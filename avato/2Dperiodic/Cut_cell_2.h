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
# ifndef CUT_CELL_2_H
# define CUT_CELL_2_H

# include <iostream>
# include "tools.h"

template<class K, class Gt>
class Cut_cell_2 // <<<
{
  public:
  typedef typename Gt::FT                   FT;
  typedef typename Gt::Point_2              Point;
  typedef typename Gt::Vector_2             Vector;
  typename K::Compute_scalar_product_2 csp;
  typename Gt::Construct_perpendicular_vector_2 cpv;

  // internal variables:
  Point E; // perpendicular drop of A on the edge

  // input:
  const Point V, // the Voronoi vertex in question
              other, // the other Vornoi vertex of the edge
              A, // sphere center (Delaunay vertex)
              another; // another Voronoi vertex of the Voronoi facet.
  const FT r2, r;

  // output:
  int sgn; // the following are without the sign:
  double vol_total;
  double vol_cavity;
  double surf;

  Cut_cell_2(const Point& V_, const Point& other_, const Point& A_, const Point& another_, const FT& r2_) // <<<
  : V(V_), other(other_), A(A_), another(another_), r2(r2_), r(CGAL::sqrt(r2_))
  {
    Vector edge(V, other);
    FT edge_lensq = squared_distance(V, other);
    Vector perpen = cpv(edge, CGAL::COUNTERCLOCKWISE);

    // calculate the perpendicular drop of A on the edge:
    if (edge_lensq > 0) {
      FT proj = csp(edge, Vector(V, A)) / edge_lensq;
      E = Point(V.x() + proj*edge.x(),
                V.y() + proj*edge.y());
    } else {
      E = V;
    }

    // Facet sign:
    // Check whether the Delaunay vertex and another Voronoi vertex are on
    // the same side of the edge: sign=+1 for same side
    // The double sf carries the sign
    FT Sf = (csp(perpen, Vector(V, A)) *
             csp(perpen, Vector(V, another)));

    // Edge sign:
    // Check whether other and the perpendicular drop E are on the same side of me(V):
    // sign=+1 for same side
    FT Se = (csp(edge, Vector(V, E)) *
             csp(edge, Vector(V, other)));

    //std::cout << signum(Sf) << " " << signum(Se)
    //          << "   " << V
    //          << "   " << other
    //          << "   " << A
    //          << "   " << another
    //          << std::endl;
    //assert(Sf >= 0); // XXX we work with good radius ratios so far ... but the code works for all others as well :-)
    sgn = CGAL::sign(Sf*Se);

    // Measure the cavity volume/surface using A,E,V,r2:
    // without applying the sign!
    // see notes p.3239
    FT x2 = CGAL::squared_distance(A, E);
    FT y2 = CGAL::squared_distance(E, V);
    FT r2_m_x2 = r2 - x2;
    //assert(std::isfinite(r2_m_x2));

    FT x = CGAL::sqrt(x2);
    FT y = CGAL::sqrt(y2);
    vol_total = CGAL::to_double(x*y/2); // total area of subcomplex

    if (r2_m_x2 >= y2) // the subcell is filled: case 3
    {
      vol_cavity = surf = 0.0;
      return;
    }
    else if (r2_m_x2 <= 0) // only a cone: case 1
    {
      double phi = (y < x ? atan(CGAL::to_double(y/x)) : M_PI_2 - atan(CGAL::to_double(x/y)));
      assert(std::isfinite(phi));
      vol_cavity = vol_total - phi*CGAL::to_double(r2/2);
      surf = CGAL::to_double(r) * phi;
    }
    else
    {
      FT yp = CGAL::sqrt(r2_m_x2); // XXX guaranteed to be real
      // we know that x^2 < r^2 but also x/r < 1?
      FT x_r = CGAL::min(FT(1), x/r);
      // avoid the M_PI/2 to gain precision:
      double phi_m_alpha = (y < x ? atan(CGAL::to_double(y/x)) - acos(CGAL::to_double(x_r))
                                  : asin(CGAL::to_double(x_r)) - atan(CGAL::to_double(x/y)));
      //assert(std::isfinite(yp));
      assert(std::isfinite(phi_m_alpha));
      vol_cavity = vol_total - (CGAL::to_double(x*yp) + CGAL::to_double(r2)*phi_m_alpha)/2;
      surf = CGAL::to_double(r) * phi_m_alpha;
    }

    // some tests
    if (vol_cavity < 0) {
        if (vol_cavity < -1.0e-6*vol_total and vol_cavity < -1.0e-12) {
          std::cerr << "WARNING: negative area vol_cavity: "
                    << vol_cavity << " "
                    << (r2 > x2)
                    << std::endl;
        }
        vol_cavity = 0;
    }
    if (surf < 0) {
        std::cerr << "WARNING: negative outline surf: "
                  << surf << " "
                  << (r2 > x2)
                  << std::endl;
        surf = 0;
    }
    if (vol_cavity > vol_total) {
        std::cerr << "WARNING: vol_cavity larger than vol_total: "
                  << vol_cavity << " " << vol_total << " "
                  << (r2 > x2)
                  << std::endl;
        vol_cavity = vol_total;
    }
  } // >>>
}; // >>>


# endif // CUT_CELL_2_H
// vim:foldmethod=marker:foldmarker=<<<,>>>
