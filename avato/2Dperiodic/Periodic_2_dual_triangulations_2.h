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

// Pair of dual tesselations in 2D: periodic regular triangulations
// - allow non-square periodic domains
// - generate Voronoi structure as additional info in the Delaunay triangulation
// - can generate several triangulations with a particle omitted and the other radii extended
// - identifies and measures cavities

#ifndef PERIODIC_2_DUAL_TRIANGULATIONS_2_H
#define PERIODIC_2_DUAL_TRIANGULATIONS_2_H

#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <bitset>
#include <algorithm>

#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Union_find.h>

# ifndef MAXLINELEN
#   define MAXLINELEN 65535
# endif
# define NDIM 2
# define SN 3
# define SN2 3 /*   number of pairs (SN * (SN-1)) // 2   */
# include "Periodic_triangulation_ds_circulators_2.h"
# include "Periodic_triangulation_cover_iterators_2.h"
# include "read_dumpfile.h"
# include "Cut_cell_2.h"
# include "tools.h"

# define HASH_WITH_IDS 1

template<class Gt, class W>
class Vor_vertex_2 // <<<
{
  public:
  typedef typename Gt::Point_2   Point;
  typedef typename Gt::Offset_2  Offset;
  Point point; // a point in the fundamental domain
  Offset offset; // offset of the point to put it at the right position,
  // relative to its dual Delaunay face.
  // NOTICE that the offsets are not limited to the cover. They can be negative, or even large...
  W weight; // weight of the Voronoi vertex
  bool is_covered; // is the Voronoi vertex covered by Spheres?
  bool dual_is_canonical; // is the dual (Delaunay face) unique?
  //int cluster_id;
  Vor_vertex_2()
  : point(0,0),
    offset(0,0),
    weight(0),
    is_covered(false),
    dual_is_canonical(false)
    //cluster_id(-1)
  {};
}; // >>>

template<class W>
class Vor_face_2 // <<<
{
  public:
  W weight; // weight of dual Delaunay vertex
  size_t unique_id; // unique ID of the Delaunay vertex (must be size_t for building unique_hash)
  //bool dual_is_fundamental;
  Vor_face_2() : weight(0), unique_id(-1) {};
}; // >>>


// REMARKS:
// * Limitation: If too dilute, the 1x1 sheet may not suffice. The whole
//   construction may fail in this case.
// * CGAL::Periodic_2_triangulation_2 needs square domain
//   and keeps private things such as
//      _too_long_edges_counter
//   Their test on _too_long_edges is stronger than our (expensive) cycle-2 test,
//   such that we cannot achieve a tri.is_valid() -- not even if we modify its private
//   data.
// * CGAL::Periodic_2_triangulation_2::get_neighbor_offset
//   does not return an offset (!), but a "wrap" in {0,1}
//   --> must use tri.combine_offsets(Offset(), tri.get_neighbor_offset(f, i));
// * CGAL::Periodic_2_triangulation_2
//   has no locate function which returns the Shift at the same time.
//   They implemented it, but make the Shift const, and protect the whole function.
//   I patched the header file to add this functionality.
//

template <class K>
class Periodic_2_dual_triangulations_2
: public CGAL::Triangulation_cw_ccw_2
{
  typedef Periodic_2_dual_triangulations_2<K>                   Self;

public:
  // typedefs <<<
  typedef CGAL::Regular_triangulation_euclidean_traits_2<K>     GtWght;
  typedef typename GtWght::Weight                               Weight;

  // we need a periodic (weighted Delaunay) triangulation which stores ts duals as well:
  typedef CGAL::Periodic_2_triangulation_traits_2<K>            GtPeri;
  typedef CGAL::Periodic_2_triangulation_vertex_base_2<GtPeri>  VbPeri;
  typedef CGAL::Periodic_2_triangulation_face_base_2<GtPeri>    FbPeri;
  typedef Vor_vertex_2<GtPeri, Weight>                          Vor_vertex;
  typedef Vor_face_2<Weight>                                    Vor_face;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vor_face, GtPeri, VbPeri> Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<Vor_vertex, GtPeri, FbPeri> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>          Tds;

  typedef CGAL::Periodic_2_triangulation_2<GtPeri, Tds>         Tri;

  // we also need to know about nonperiodic regular triangulations:
  typedef CGAL::Regular_triangulation_vertex_base_2<GtWght>     VbWght;
  typedef CGAL::Regular_triangulation_face_base_2<GtWght>       FbWght;
  typedef CGAL::Triangulation_data_structure_2<VbWght, FbWght>  TdsReg;
  typedef CGAL::Regular_triangulation_2<GtWght, TdsReg>         Reg;

  typedef typename GtPeri::Periodic_2_offset_2  Offset;
  typedef Offset                                Shift; // the difference between two Offsets
  typedef typename GtPeri::Iso_rectangle_2      Iso_rectangle;
  typedef typename GtPeri::FT                   FT;
  typedef typename GtPeri::Point_2              Point;
  typedef typename GtPeri::Vector_2             Vector;
  typedef typename GtWght::Weighted_point_2     Weighted_point;
  typedef typename K::Compute_scalar_product_2 Compute_scalar_product;
  typedef typename GtPeri::Construct_perpendicular_vector_2 Construct_perpendicular_vector;

  typedef typename Tri::Vertex_handle                   Del_Vertex_handle;
  typedef typename Tri::Vertex_iterator                 Del_Vertex_iterator;
  typedef typename Tri::Unique_vertex_iterator          Del_Unique_vertex_iterator;
  typedef typename Tri::Face_handle                     Del_Face_handle;
  typedef typename Tri::Face_iterator                   Del_Face_iterator;
  typedef Periodic_triangulation_ds_face_circulator_2<Tds> Del_My_Face_circulator;
  typedef Periodic_face_cover_iterator_2<Tri>           Del_Face_cover_iterator;

  typedef CGAL::Handle_hash_function            Hhash;
  //pedef Handle_hash_function                  Hhash;
  typedef Handle_Offset_2_hash_function<Hhash>  HOhash;
  typedef CGAL::Union_find<Del_Face_handle>                     Clusters;
  typedef typename Clusters::handle                             Cluster;
  typedef std::unordered_map<Del_Face_handle, Cluster, Hhash>   DF_TO_CLS;
  typedef Cut_cell_2<K, GtPeri>                                 Cut_cell;
  typedef Cavity_measure_base<FT>                               Cavity_measure;
  typedef std::unordered_map<Cluster, Cavity_measure, Ptr_hash_function> Cavity_measures;
  // >>>

private:
  typedef typename Reg::Vertex_handle                   Reg_Vertex_handle;
  typedef typename Reg::Face_handle                     Reg_Face_handle;
  typedef std::unordered_map<Reg_Vertex_handle, std::pair<Reg_Vertex_handle, Offset>, Hhash>    REGV_TO_REGV_OFF;
  typedef std::unordered_map<Reg_Vertex_handle, std::pair<Del_Vertex_handle, Offset>, Hhash>    REGV_TO_TRIV_OFF;
  typedef std::unordered_map<std::pair<Del_Vertex_handle, Offset>, Del_Vertex_handle, HOhash>   TRIV_OFF_TO_TRIV;
  typedef std::unordered_map<Del_Face_handle, Reg_Face_handle, Hhash>                           TRIF_TO_REGF;
  typedef std::unordered_map<Reg_Face_handle, Del_Face_handle, Hhash>                           REGF_TO_TRIF;
  typedef std::unordered_map<Reg_Face_handle, std::pair<Del_Face_handle, int>, Hhash>           REGF_TO_TRIF_CYCLE;
  typedef std::unordered_map<Del_Face_handle, std::pair<Del_Face_handle, Shift>, Hhash>         TRIF_TO_TRIF_SH;
  typedef std::unordered_map<std::pair<Del_Face_handle, Shift>, Del_Face_handle, HOhash>        TRIF_SH_TO_TRIF;
  typedef std::map<size_t, Del_Face_handle>                                                     HASH_TO_TRIF;
  typedef typename REGV_TO_REGV_OFF::const_iterator     REGV_TO_REGV_OFF_cit;
  typedef typename REGV_TO_TRIV_OFF::const_iterator     REGV_TO_TRIV_OFF_cit;
  typedef typename TRIV_OFF_TO_TRIV::const_iterator     TRIV_OFF_TO_TRIV_cit;
  typedef typename TRIF_TO_REGF::const_iterator         TRIF_TO_REGF_cit;
  typedef typename REGF_TO_TRIF::const_iterator         REGF_TO_TRIF_cit;
  typedef typename REGF_TO_TRIF::iterator               REGF_TO_TRIF_it;
  typedef typename REGF_TO_TRIF_CYCLE::const_iterator   REGF_TO_TRIF_CYCLE_cit;
  typedef typename REGF_TO_TRIF_CYCLE::iterator         REGF_TO_TRIF_CYCLE_it;
  typedef typename TRIF_TO_TRIF_SH::const_iterator      TRIF_TO_TRIF_SH_cit;
  typedef typename TRIF_SH_TO_TRIF::const_iterator      TRIF_SH_TO_TRIF_cit;
  typedef typename HASH_TO_TRIF::const_iterator         HASH_TO_TRIF_cit;

  static const size_t SIZE_T_NONE = SIZE_MAX-10;

  Iso_rectangle _domain;
  const FT PREC;
  std::vector<Point> _orig_points;
  std::vector<double> _orig_radii;
  double _min_radius, _max_radius;
  bool _one_sheet_is_always_sufficient = false;
  bool _nine_sheets_are_always_sufficient = false;
  FT _halfway_squared = -1.0;

public:

// Default constructor
  Periodic_2_dual_triangulations_2(const Iso_rectangle& domain = Iso_rectangle(0, 0, 1, 1), const FT& prec=1.0e-10)
    : _domain(domain), PREC(prec)
  {}

// user interface
  const Iso_rectangle& domain() const {return _domain;}
  const std::vector<Point>& orig_points() const {return _orig_points;}
  const std::vector<double>& orig_radii() const {return _orig_radii;}
  void overwrite_radii(const std::vector<double>& radii) {_orig_radii = radii;}
  bool read_points(const char* filename, bool keep_domain=false) // <<<
  {
    _orig_points.clear();
    _orig_radii.clear();
    Iso_rectangle domain;
    if (not read_dumpfile<Iso_rectangle, Point>(domain, _orig_points, _orig_radii, filename)) {
      return false;
    }
    assert(_orig_points.size() == _orig_radii.size());
    if (_orig_radii.size() == 0) {
      std::cerr << "ERROR: No points in \"" << filename << "\"" << std::endl;
      exit(EXIT_FAILURE);
    }

    if (keep_domain) {
      assert(CGAL::abs(domain.xmin() - _domain.xmin()) < PREC);
      assert(CGAL::abs(domain.xmax() - _domain.xmax()) < PREC);
      assert(CGAL::abs(domain.ymin() - _domain.ymin()) < PREC);
      assert(CGAL::abs(domain.ymax() - _domain.ymax()) < PREC);
    } else {
      _domain = domain;
    }

    // half the minimal distance in the fundamental domain
    _halfway_squared = FT(0.25) * std::min(
      (_domain.xmax()-_domain.xmin()) * (_domain.xmax()-_domain.xmin()),
      (_domain.ymax()-_domain.ymin()) * (_domain.ymax()-_domain.ymin()));

    // wrap _orig_points into the domain:
    for (size_t n=0, N=_orig_points.size(); n<N; n++) {
      _orig_points[n] = wrapped(_orig_points[n]);
    }

    // min_radius and max_radius
    typedef const std::vector<double>::const_iterator it_type;
    std::pair<it_type, it_type> mm = std::minmax_element(_orig_radii.cbegin(), _orig_radii.cend());
    _min_radius = *mm.first;
    _max_radius = *mm.second;
    assert(0 < _min_radius and _min_radius <= _max_radius);

    // see notes p.3294 and Caroli+Teillaud(2009):
    // we can guarantee that all point insertions keep the
    // triangulation simplicial, only starting from
    //   5x5 sheets in 3D
    //   4x4 sheets in 2D
    // These values are increased further if the domain
    // is not square, and if the weights are not all equal.
    _one_sheet_is_always_sufficient = false;
    _nine_sheets_are_always_sufficient = false;

    return true;
  } // >>>
  void update_orig_point(const size_t n, const Point& pt) // <<<
  {
    assert(_orig_points.size() > 0);
    _orig_points[n] = wrapped(pt, 1.0e-14);
  } // >>>
  Tri del_from_nonperiodic_regular(const size_t exclude_n=SIZE_T_NONE, const double extra_radius = -1.0) const // <<<
  {
    assert(_orig_points.size() == _orig_radii.size());
    assert(typeid(typename Reg::Bare_point) == typeid(Point));

    // mapping to fundamental set of vertices:
    REGV_TO_REGV_OFF regv_to_fundregv_off;

    // Create several copies of the regular triangulation
    Reg reg;
    for (size_t n=0, N=_orig_points.size(); n<N; n++) {
      // skip one particle
      if (n == exclude_n)
        continue;
      assert(_domain.has_on_bounded_side(_orig_points[n]));
      // enlarge the others
      Weight wght;
      if (extra_radius > 0) {
        wght = (_orig_radii[n]+extra_radius) * (_orig_radii[n]+extra_radius);
      } else {
        wght = _orig_radii[n] * _orig_radii[n];
      }
      // insert point into regular triangulation
      Reg_Vertex_handle fundv = reg.insert(Weighted_point(_orig_points[n], wght));
      regv_to_fundregv_off[fundv] = std::make_pair(fundv, Offset(0,0));

      // According to Dolbilin+Huson(1996) 3x3 copies are always sufficient
      // to see a periodic triangulation (which is not necessarily simplicial
      // in the sense of Caroli+Teillaud(2009))
      // see notes p.3298 for weighted Delaunay triangulation
      for (int i=-1; i<2; i++) {
        for (int j=-1; j<2; j++) {
          if (i == 0 and j == 0) continue;
          Offset off(i,j);
          Reg_Vertex_handle v = reg.insert(Weighted_point(unwrapped(_orig_points[n], off), wght));
          regv_to_fundregv_off[v] = std::make_pair(fundv, off);
        }
      }
    }
    assert((exclude_n != SIZE_T_NONE) or (reg.number_of_vertices() + reg.number_of_hidden_vertices() == 9*_orig_points.size()));
    assert(reg.is_valid());

    // draw:
    //std::cout << "Original points: " << _orig_points.size() << std::endl;
    //std::cout << "Regular vertices, hidden vertices: " << reg.number_of_vertices() << ", " << reg.number_of_hidden_vertices() << std::endl;
    //reg_save("regular.nptri", reg);

    // create our periodic weighted Delaunay triangulation from the regular one
    // CAUTION! the domain is non-square, we do the necessary checks ourselves.
    // From a CGAL point of view, tri is *invalid* in the following sense:
    // * it does not have correct _too_long_edges
    Tri tri;
    if (reg.number_of_vertices() == 0) return tri;
    tri.set_domain(_domain); // implies clear()
    tri.tds().set_dimension(2); // must be done explicitly
    // we do the main work in a subroutine:
    del_create_1_sheeted(tri, reg, regv_to_fundregv_off);
    //del_comply_internals(tri);
    # ifndef NDEBUG
    if (not _one_sheet_is_always_sufficient)
      if (contains_2_cycle(tri))
      {
        # if 0
        std::cerr << "WARNING: Detected 2-cycle. Going on with 1 sheet, fingers crossed..." << std::endl;
        # else
        std::cerr << "INFO: Need 3x3 covering for periodic weighted Delaunay triangulation." << std::endl;
        tri.convert_to_9_sheeted_covering();
        if (not _nine_sheets_are_always_sufficient)
          if (contains_2_cycle(tri))
          {
            std::cerr << "ERROR: Need more than 3x3 sheets for periodic weighted Delaunay triangulation." << std::endl;
            // TODO: raise exception instead of exit
            exit(EXIT_FAILURE);
          }
        # endif
      }
    //assert(tri.is_valid()); // Cannot be asserted here,
    // because above we check for 2-cycles, not for tri._too_long_edges.empty().
    # endif
    # if HASH_WITH_IDS
    if (tri.number_of_vertices() >= (1<<16)) {
      std::cout << "Too many vertices for del_unique_face_hash. Try CGAL handle hash function." << std::endl;
      exit(EXIT_FAILURE);
    }
    # endif
    return tri;
  } // >>>
  void vor_vertices_create(Tri& tri, __attribute__((unused)) const FT& prec=1.0e-10) const // <<<
  // create Voronoi structure in the info() of Delaunay faces
  // sets:
  //      vv.weight
  //      vv.point
  //      vv.offset
  //      vv.is_covered
  {
    if (tri.empty()) return;

    Weight max_radius_squared = del_largest_weight(tri);
    const Reg reg; // temporary regular to construct circumcenter
    std::map<size_t, Point> images;
    size_t n_canonicals = 0;

    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it; // cast
      Vor_vertex& vv = df->info();
      // vv.dual_is_canonical was set at construction of the Delaunay Facets
      if (vv.dual_is_canonical) n_canonicals++;

      // set point and offset:
      vv.point = wrapped(vv.offset, weighted_circumcenter(vv.weight, df, tri, reg));
      assert(tri.domain().has_on_bounded_side(vv.point));

      // reset point: avoid different copies of the same point
      size_t hash = del_unique_face_hash(df, tri);
      if (images.find(hash) == images.end()) {
        images[hash] = vv.point;
      } else {
        assert(squared_distance(vv.point, images[hash]) < prec);
        vv.point = images[hash];
      }

      // is_covered:
      // Search (with zero-offset) the Delaunay face which contains the point.
      // This can be different from df
      typename Tri::Locate_type lt;
      int li;
      Offset loc_off;
      Del_Face_handle loc_df = tri.locate(vv.point, loc_off, lt, li, df);
      assert(loc_df != Del_Face_handle());
      vv.is_covered = vor_vertex_is_covered_by_spheres(vv, loc_df, -loc_off, tri, max_radius_squared);
    }
    // make sure we found the right number of canonical faces
    assert(n_canonicals == images.size() and images.size() == tri.number_of_faces());

  } // >>>
  bool vor_faces_contain_duals(const Tri& tri) const // <<<
  {
    for (Del_Unique_vertex_iterator dv_it = tri.unique_vertices_begin();
         dv_it != tri.unique_vertices_end();
         ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      assert(tri.get_offset(dv).is_zero());
      if (not vor_face_contains_point(dv->point(), dv, tri.get_offset(dv), tri))
        return false;
    }
    return true;
  } // >>>
  std::pair<Del_Vertex_handle, Offset> vor_locate(const Point& point, const Tri& tri) const // <<<
  // Returns the Voronoi face (via its dual) which contains the point p. The
  // second argument is the offset which must be applied instead of
  // tri.get_offset(dv) to the Voronoi face (Delaunay vertex) in order to
  // contain the point.
  //
  // Algorithm:
  // (We can never exclude that a point on an edge is rejected by both facets:
  //  Predicates are dependent on the reconstructed periodic points, thus they are never exact!)
  //  1. First, find the closest Voronoi vertex:
  //     1a Do a march on the Voronoi edges until local distance minimum
  //     1b Check the three incident Voronoi facets for a closer vertex.
  //        If found one, continue with 1a
  //  2. Now, the point *must* be in one of the three incident Voronoi facets.
  //     Use a must-have-answer function.
  {

    Shift df_shift;
    Del_Face_handle df = tri.faces_begin();
    while(true) {
      // 1a
      vor_locate_edge_march(df, df_shift, point, tri);
      // 1b
      if (not vor_locate_incident_faces(df, df_shift, point, tri)) break;
    }
    // the closest periodic copy of a Voronoi vertex is now
    //   unwrapped(dc->info().point, dc->info().offset + dc_shift)

    // 2
    int take_dv_i = vor_vertex_sector_of_point(df, df_shift, point, tri);
    return std::make_pair(df->vertex(take_dv_i),
                          tri.get_offset(df, take_dv_i) + df_shift);
  } // >>>
  bool vor_locate_edge_march(Del_Face_handle& df, Shift& df_shift, const Point& point, const Tri& tri) const // <<<
  // The Voronoi vertex (df) follows the incident Voronoi edges until it finds a local distance minimum to point.
  // NOTE: To avoid ambiguities, there is only *one* copy of point, but df can move over all images by means of df_shift.
  // Return value:
  //   true if we moved
  //   false if df,df_shift is already the local minimum
  // Output:
  //   df the new closest Voronoi vertex
  //   df_shift it additional shift to be put at the correct position
  {
    const Vor_vertex& vv = df->info();
    FT distsq = squared_distance(point, unwrapped(vv.point, vv.offset + df_shift));

    int came_from = -1;

    while (true) {

      bool found_closer = false;
      Del_Face_handle closer_df = df;
      Shift closer_df_shift = df_shift;

      // find the closest of the neighbors
      for (int i=0; i<3; i++) {
        if (i == came_from) continue;
        const Vor_vertex& nbr_vv = df->neighbor(i)->info();
        Shift nbr_df_shift = tri.combine_offsets(df_shift, -tri.get_neighbor_offset(df, i));
        FT nbr_distsq = squared_distance(point, unwrapped(nbr_vv.point, nbr_vv.offset + nbr_df_shift));
        if (nbr_distsq <= distsq) {
          // We can use <= because we never go back the path we came from
          // It may happen that two VorVertices have equal distance from the point,
          // but only one of them knows about closer points...
          found_closer = true;
          closer_df = df->neighbor(i);
          closer_df_shift = nbr_df_shift;
          distsq = nbr_distsq;
        }
      }

      if (found_closer) {
        came_from = closer_df->index(df);
        df = closer_df;
        df_shift = closer_df_shift;
        continue;
      }
      return (came_from != -1);
    }
    assert(false); // cannot come here
    return true;
  } // >>>
  bool vor_locate_incident_faces(Del_Face_handle& df, Shift& df_shift, const Point& point, const Tri& tri) const // <<<
  // The Voronoi vertex (df) is set to the VorVertex which is closest (<) to
  // point, choosing among the VorVertices of all incident VorFaces.
  //
  // NOTE: To avoid ambiguities, there is only *one* copy of point, but df can
  // move over all images by means of df_shift.
  //
  // Return value:
  //   true if we moved
  //   false if df,df_shift is already the local minimum
  // Output:
  //   df the new closest Voronoi vertex
  //   df_shift it additional shift to be put at the correct position
  {
    Del_Face_handle closest_df = df;
    Shift closest_df_shift = df_shift;
    FT closest_distsq;
    {
      const Vor_vertex& vv = df->info();
      closest_distsq = squared_distance(point, unwrapped(vv.point, vv.offset + df_shift));
    }

    // Loop over the vertices of df (== the incident VorFaces of VorVertex(df))
    for (int i=0; i<3; i++) {
      Del_Vertex_handle dv = df->vertex(i);

      // Circulate around the vertices of the VorFace(dv),
      // starting at df
      const Del_My_Face_circulator df_it(dv, df, i);
      Del_My_Face_circulator nbr_df_it(df_it);
      // TODO: avoid df and the following point (recurrent calculation of distances)
      do {
        Del_Face_handle nbr_df = nbr_df_it;
        const Vor_vertex& nbr_vv = nbr_df->info();
        // df and nbr_df have a common vertex. The new face shift is enlarged by the shift of this vertex:
        Shift nbr_df_shift = df_shift - tri.get_offset(nbr_df, nbr_df_it.index()) + tri.get_offset(df, i);
        FT nbr_distsq = squared_distance(point, unwrapped(nbr_vv.point, nbr_vv.offset + nbr_df_shift));
        if (nbr_distsq < closest_distsq) {
          // Use < here to avoid loops between vor_locate_edge_march and vor_locate_incident_faces
          closest_df = nbr_df;
          closest_df_shift = nbr_df_shift;
          closest_distsq = nbr_distsq;
        }
      } while (++nbr_df_it != df_it);

    }

    bool did_move = (closest_df != df or closest_df_shift != df_shift);
    df = closest_df;
    df_shift = closest_df_shift;
    return did_move;
  } // >>>
  void vor_find_all_clusters(Clusters& clusters, DF_TO_CLS& df_to_cls, const Tri& tri) const // <<<
  // Finds connectivity clusters of Voronoi vertices.
  // Returns the Union_find object and
  // the mapping Delaunay faces->clusters in the Union_find
  // On return:
  //   * Clusters contains only non-covered Voronoi vertices which have a canonical dual
  //   * df_to_cls is a mapping from _all_ Delaunay faces to the cluster of their dual
  //   * neither clusters nor df_to_cls contain any covered Vor_vertex
  //
  // For later use:
  // The iterator of Clusters runs over all elements.
  // We can identify the class by ptr() of the iterator.
  // XXX We do this on the fundamental domain only!
  {
    clusters.clear();
    df_to_cls.clear();
    if (tri.empty()) return;
    if (not tri.is_1_cover())
      std::cerr << "INFO: Triangulation is not 1-sheeted. Extracting clusters only from the fundamental domain..." << std::endl;

    // we need to collect the canonical faces
    std::map<size_t, Del_Face_handle> canon;

    // XXX we have already set the "is_covered" property on all Voronoi vertices
    // 1. Collect all non-covered Voronoi vertices
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      const Del_Face_handle df = df_it; // cast operator
      const Vor_vertex& vv = df->info();
      if (vv.dual_is_canonical) {
        size_t hash = del_unique_face_hash(df, tri);
        assert(canon.find(hash) == canon.end());
        canon.insert(std::make_pair(hash, df));
        if (not vv.is_covered)
          df_to_cls[df] = clusters.make_set(df);
      }
    }

    // 2. identify clusters by checking whether the edges are "covered"
    // XXX in this loop, every Voronoi edge is checked twice
    for (typename Clusters::iterator cl_it = clusters.begin();
         cl_it != clusters.end();
         ++cl_it)
    {
      const Del_Face_handle df = *cl_it;
      const Vor_vertex& vv = df->info();
      const Point vp = unwrapped(vv.point, vv.offset); // point with respect to df

      // Circulate over Voronoi neighbor vertices
      // by means over the Delaunay neighbor 
      for (int i=0; i<3; i++) {
        Del_Face_handle df_nbr = df->neighbor(i);
        // df_nbr must consider the shift from df to df->neighbor(i)
        Shift nbr_shift = tri.get_offset(df, ccw(i)) - tri.get_offset(df_nbr, cw(df_nbr->index(df)));

        if (not df_nbr->info().dual_is_canonical) {
          // we need a different nbr (and different vv_nbr below)
          // which means we must also correct its vv_nbr.offset
          size_t hash = del_unique_face_hash(df_nbr, tri);
          assert(canon.find(hash) != canon.end());
          Del_Face_handle canon_nbr = canon[hash];
          assert(tri.get_original_vertex(df_nbr->vertex(0)) ==
                 tri.get_original_vertex(canon_nbr->vertex(0)));
          nbr_shift += tri.get_offset(df_nbr, 0) - tri.get_offset(canon_nbr, 0);
          df_nbr = canon_nbr;
        }

        const Vor_vertex& vv_nbr = df_nbr->info();
        Del_Vertex_handle dvA = df->vertex(ccw(i));
        // dvA and dvB are both vertices with respect to the face df:
        Point dpA = unwrapped(dvA->point(), tri.get_offset(df, ccw(i))),
           vp_nbr = unwrapped(vv_nbr.point, vv_nbr.offset + nbr_shift);
        bool resA = vor_edge_is_split_by_sphere(vp, vp_nbr, dpA, dvA->info().weight);

        # ifndef NDEBUG
        // we have two spheres at our choice -- make sure they give the same answer.
        Del_Vertex_handle dvB = df->vertex(cw(i));
        Point dpB = unwrapped(dvB->point(), tri.get_offset(df, cw(i)));
        bool resB = vor_edge_is_split_by_sphere(vp, vp_nbr, dpB, dvB->info().weight);
        if (resA != resB) {
          std::cerr
            << vp << " \t" << vp_nbr << std::endl
            << dpA << " \t" << dvA->info().weight << std::endl
            << dpB << " \t" << dvB->info().weight << std::endl;
        }
        assert(resA == resB);
        # endif

        if (not resA) {
          assert(clusters.find(df_to_cls[df]) != clusters.end());
          assert(clusters.find(df_to_cls[df_nbr]) != clusters.end());
          clusters.unify_sets(clusters.find(df_to_cls[df]), clusters.find(df_to_cls[df_nbr]));
          assert(clusters.find(df_to_cls[df]) == clusters.find(df_to_cls[df_nbr]));
          //assert(df_to_cls[df] == df_to_cls[df_nbr]); is not the case!
        }
      }
    }

    // For later use, we extend df_to_cls to all faces
    // and we pre-search the clusters
    if (not tri.is_1_cover())
    {
      for (Del_Face_iterator df_it = tri.faces_begin();
           df_it != tri.faces_end();
           ++df_it)
      {
        Del_Face_handle df = df_it; // cast operator
        if (df->info().is_covered)
          continue;

        Cluster cl = df_to_cls[df];
        if (not df->info().dual_is_canonical) {
          size_t hash = del_unique_face_hash(df, tri);
          assert(canon.find(hash) != canon.end());
          cl = df_to_cls[canon[hash]];
        }
        df_to_cls[df] = clusters.find(cl);
      }
    }

    # if 0
    // set the cluster_id in the Voronoi vertices
    std::unordered_map<Cluster, size_t, Ptr_hash_function> cl_to_int;
    int n = 0;
    for (typename Clusters::iterator cl_it = clusters.begin();
         cl_it != clusters.end();
         ++cl_it)
    {
      Del_Face_handle df = *cl_it;
      assert(not df->info().is_covered);
      Cluster cl = clusters.find(cl_it);
      if (cl_to_int.find(cl) == cl_to_int.end())
        cl_to_int[cl] = n++;
      df->info().cluster_id = cl_to_int[cl];
    }
    # endif

    # if 0
    std::cerr << "Cluster membership:" << std::endl;
    for (typename Clusters::iterator cl_it = clusters.begin();
         cl_it != clusters.end();
         ++cl_it)
    {
      Del_Face_handle df = *cl_it;
      std::cerr
        << " " << clusters.find(cl_it).ptr() // important! cl_itr.ptr() is not the same
        << " \t" << df->info().point
        << std::endl;
    }
    # endif

  } // >>>
  Cluster vor_locate_in_clusters(const Point& point, Clusters& clusters, DF_TO_CLS& df_to_cls, const Tri& tri, const FT& prec=1.0e-10) const // <<<
  // Returns the cluster which contains the point (or Cluster(NULL) if point is in no cluster)
  // prec is a precision parameter of dimension length^2
  //
  // Algorithm:
  // 1. First find the Voronoi facet containing the point.
  // 2. Then identify the corresponding Voronoi vertex in that VorFacet -- and its cluster in which the point is located.
  // TODO: make the sphere comparison a predicate.
  // TODO: how is const to be treated in Union_find?
  {
    assert(clusters.size() == vor_vertices_number_of_not_covered(tri));
    assert(df_to_cls.size() == clusters.size() * tri.number_of_sheets()[0]*tri.number_of_sheets()[1]);
    assert(prec >= 0);

    // This is the Voronoi facet (Delaunay vertex dual) and its explicit offset (may be shifted)
    std::pair<Del_Vertex_handle, Offset> loc = vor_locate(point, tri);
    Del_Vertex_handle dv = loc.first;
    const Offset& dv_off = loc.second;
    Point dp = unwrapped(dv->point(), dv_off);

    // If the point is covered by the sphere, return immediately
    // with the missed parameter set
    FT rsq_missed = dv->info().weight - squared_distance(dp, point);
    if (rsq_missed > prec) return Cluster(NULL);

    // Identify the VorVertex in the VorFacet
    //Del_My_Face_circulator df_it = vor_face_sector_containing_point(dv, dv_off, point, tri, prec); // old version
    Del_My_Face_circulator df_it = del_vertex_sector_containing_point(dv, dv_off, point, tri, prec);
    Del_Face_handle df = df_it; // cast
    assert(df->vertex(df_it.index()) == dv);

    if (not df->info().is_covered) {
      assert(df_to_cls.find(df) != df_to_cls.end());
      assert(clusters.find(df_to_cls[df]) != clusters.end());
      return clusters.find(df_to_cls[df]);
    }
    assert(df_to_cls.find(df) == df_to_cls.end());
    return Cluster(NULL);
  } // >>>
  Cavity_measures vor_cavity_measures(const Tri& tri, Clusters& clusters, Cluster only_cluster=Cluster()) const // <<<
  {
    Cavity_measures result;
    bool only_one = (only_cluster != Cluster());

    // loop over cluster vertices and sum up contributions
    for (typename Clusters::iterator cl_it = clusters.begin();
         cl_it != clusters.end();
         ++cl_it)
    {
      Cluster cl = clusters.find(cl_it).ptr();
      if (only_one)
        if (cl != only_cluster)
          continue;

      Cavity_measure& cm = result[cl];
      Del_Face_handle df = *cl_it;
      const Vor_vertex& vv = df->info();
      Point V = unwrapped(vv.point, vv.offset);

      // We circulate around the Voronoi vertex (Delaunay face) and collect its contributions.
      // Each Voronoi edge (Delaunay facet nbr) has two spheres to consider.
      for (int i=0; i<3; i++)
      {
        // get the other point of the edge
        const Vor_vertex& vv_other = df->neighbor(i)->info();
        Point vp_other = unwrapped(vv_other.point,
                                   tri.combine_offsets(vv_other.offset, -tri.get_neighbor_offset(df, i)));

        // consider both spheres:
        int j = ccw(i), k = cw(i);
        for (int foo=0; foo<2; foo++) // do two loops
        {
          // sphere
          Del_Vertex_handle dv = df->vertex(j);
          Point A = unwrapped(dv->point(), tri.get_offset(df, j));
          // another point of the Voronoi facet
          Del_Face_handle df_another = df->neighbor(k);
          const Vor_vertex& vv_another = df_another->info();
          Point vp_another = unwrapped(vv_another.point,
                                       tri.combine_offsets(vv_another.offset, -tri.get_neighbor_offset(df, k)));
          // add up the measure of the cavity
          Cut_cell cut(V, vp_other, A, vp_another, dv->info().weight);
          cm.add(cut.sgn*cut.vol_total,
                 cut.sgn*cut.vol_cavity,
                 cut.sgn*cut.surf,
                 cut.r);
          j = cw(i), k = ccw(i);
        } // end loop over spheres
      } // end loop over neighbors (Voronoi edges)
    }

    # if 0
    for (typename Cavity_measures::const_iterator cm_it = result.cbegin();
         cm_it != result.cend(); ++cm_it)
    {
      Cluster cl = cm_it->first;
      Cavity_measure cm = cm_it->second;
      std::cerr << clusters.find(cl).ptr()
                << " " << cm.vol_triangles
                << " " << cm.vol_cavity
                << " " << cm.surf
                << " " << cm.surf_radii
                << std::endl;
    }
    # endif
    return result;
  } // >>>

// one level higher: treat files with several triangulations
  typedef std::vector<std::pair<std::pair<size_t, FT>, Cavity_measures> > avato_t;
  avato_t treat_available_volumes_after_takeout(bool only_one_cavity, bool verbose=true, std::ostream& out=std::cerr, const FT accept_missed_cavity=1.0e-12, const FT accept_missed_cavity_max=1.0e4, __attribute__ ((unused)) bool allow_locate_to_fail=false) const // <<<
  {
    Clusters clusters;
    DF_TO_CLS df_to_cls;

    std::ios::fmtflags out_oldflags(out.flags()); // save state
    if (verbose) {
      out.precision(12);
      out << std::scientific;
    }

    avato_t result;
    for (size_t n=0; n<_orig_points.size(); n++)
    {
      if (verbose) {
        std::cerr << "INFO: Omitting particle " << n << std::endl;
        if (out.rdbuf() != std::cerr.rdbuf())
          out << "INFO: Omitting particle " << n << std::endl;
      }
      // create Delaunay and Voronoi:
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      //tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      //char filename[MAXLINELEN];
      //snprintf(filename, MAXLINELEN, "delaunay_%lu.mytri", n); del_save(filename, tri);
      //snprintf(filename, MAXLINELEN, "voronoi_%lu.mytri", n); vor_save(filename, tri);

      // identify all cavities
      vor_find_all_clusters(clusters, df_to_cls, tri);
      assert((clusters.number_of_sets() >= 1 and clusters.size() >= 1) or allow_locate_to_fail);
      if (verbose)
        std::cerr << "INFO: Clusters: " << clusters.size() << " Voronoi vertices in " << clusters.number_of_sets() << " sets." << std::endl;

      Cluster cl = Cluster(NULL);
      // find the cavity with the omitted center in it
      FT prec = accept_missed_cavity;
      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri, prec);
      while (cl == Cluster(NULL) and prec < accept_missed_cavity_max) {
        prec *= 2;
        if (verbose)
          std::cerr << "WARNING: Point " << _orig_points[n] << " not found in the only cluster. Repeating search with prec=" << prec << std::endl;
        cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri, prec);
      }
      assert(cl != Cluster(NULL) or allow_locate_to_fail);

      Cavity_measures cms;
      if (only_one_cavity) {
        if (cl != Cluster(NULL)) {
          // measure this cavity
          cms = vor_cavity_measures(tri, clusters, cl);
          assert(cms.size() == 1 or allow_locate_to_fail);
        }
      } else {
        // measure all cavities
        cms = vor_cavity_measures(tri, clusters);
      }
      result.push_back(std::make_pair(std::make_pair(n, _orig_radii[n]), cms));

      // print
      if (verbose)
      {
        for (typename Cavity_measures::const_iterator cm_it = cms.begin();
             cm_it != cms.end();
             ++cm_it)
        {
          const Cavity_measure& cm = cm_it->second;
          out << "  " << n << " OK"
              << " " << _orig_radii[n]
              << " " << cm.vol_cavity
              << " " << cm.surf
              << " " << cm.vol_triangles
              << " " << cm.surf_radii;
          if (cm_it->first == cl) out << " *";
          out << std::endl;
        }
      }
    }
    if (verbose)
      out.flags(out_oldflags); // restore flags
    return result;
  } // >>>
  typedef std::vector<std::pair<FT, Cavity_measures> > av_t;
  av_t treat_available_volumes(bool verbose=true, std::ostream& out=std::cerr, const double prec=1.0e-5) const // <<<
  {
    Clusters clusters;
    DF_TO_CLS df_to_cls;

    std::ios::fmtflags out_oldflags(out.flags()); // save state
    if (verbose) {
      out.precision(12);
      out << std::scientific;
    }

    // get radius classes
    std::vector<double> radii;
    if (prec > 0.0)
    {
      // get radius classes (scales at worst N^2)
      // the only guarantee is that representatives are at least prec apart
      // the resulting binning is not uniform.
      std::set<double> tmp(_orig_radii.cbegin(), _orig_radii.cend());
      const double precsq = prec*prec;
      for (std::set<double>::const_iterator orig = tmp.cbegin(); orig != tmp.cend(); ++orig)
      {
        bool addthis = true;
        for (std::vector<double>::const_iterator r = radii.cbegin(); r != radii.cend(); ++r)
        {
          if ((*orig - *r)*(*orig - *r) < precsq)
          {
            addthis = false;
            break;
          }
        }
        if (addthis)
          radii.push_back(*orig);
      }
      std::sort(radii.begin(), radii.end());
    }
    else
    {
      // treat all radii as their own class
      radii = _orig_radii;
    }


    av_t result;
    size_t n = 0;
    for (std::vector<double>::const_iterator r = radii.cbegin();
         r != radii.cend(); ++n, ++r)
    {
      if (verbose) {
        std::cerr << "INFO: Treating radius #" << n << ": " << *r << std::endl;
        if (out.rdbuf() != std::cerr.rdbuf())
          out << "INFO: Treating radius #" << n << ": " << *r << std::endl;
      }
      // create Delaunay and Voronoi:
      Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, *r);
      //tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);

      // identify all cavities
      vor_find_all_clusters(clusters, df_to_cls, tri);
      if (verbose)
        std::cerr << "INFO: Clusters: " << clusters.size() << " Voronoi vertices in " << clusters.number_of_sets() << " sets." << std::endl;
      // measure the cavities
      result.push_back(std::make_pair(*r, vor_cavity_measures(tri, clusters)));

      // print
      if (verbose)
      {
        Cavity_measures& cms = result[result.size()-1].second;
        for (typename Cavity_measures::const_iterator cm_it = cms.begin();
             cm_it != cms.end();
             ++cm_it)
        {
          const Cavity_measure& cm = cm_it->second;
          out << "  " << n << " OK"
              << " " << *r
              << " " << cm.vol_cavity
              << " " << cm.surf
              << " " << cm.vol_triangles
              << " " << cm.surf_radii
              << std::endl;
        }
      }
    }
    if (verbose)
      out.flags(out_oldflags); // restore flags
    return result;
  } // >>>

private:
  void del_create_1_sheeted(Tri& tri, const Reg& reg, const REGV_TO_REGV_OFF& regv_to_fundregv_off) const // <<<
  {
    typedef typename Reg::Finite_vertices_iterator        Reg_Finite_vertices_iterator;
    typedef typename Reg::Finite_faces_iterator           Reg_Finite_faces_iterator;
    REGV_TO_TRIV_OFF regv_to_triv_off;

    // VERTICES
    // vertices in the 1-sheet domain (which coincides with _domain)
    // XXX "finite" here means "not hidden"
    size_t unique_id = 0;
    for (Reg_Finite_vertices_iterator regv_fit = reg.finite_vertices_begin();
         regv_fit != reg.finite_vertices_end();
         ++regv_fit)
    {
      Reg_Vertex_handle regv = regv_fit; // cast operator
      REGV_TO_REGV_OFF_cit fundregv_it = regv_to_fundregv_off.find(regv);
      assert(fundregv_it != regv_to_fundregv_off.end());
      if (fundregv_it->second.first == regv) {
        // create one vertex
        Del_Vertex_handle triv = tri.tds().create_vertex();
        triv->set_point(regv->point().point());
        triv->info().weight = regv->point().weight();
        triv->info().unique_id = unique_id++;
        // store the mapping regv -> triv + offset
        regv_to_triv_off[regv] = std::make_pair(triv, Offset());
      }
    }
    // remaining vertices
    for (Reg_Finite_vertices_iterator regv_fit = reg.finite_vertices_begin();
         regv_fit != reg.finite_vertices_end();
         ++regv_fit)
    {
      Reg_Vertex_handle regv = regv_fit; // cast operator
      REGV_TO_TRIV_OFF_cit triv_it = regv_to_triv_off.find(regv);
      if (triv_it == regv_to_triv_off.end())
      {
        // some of the outer images may have hidden fundamental images:
        REGV_TO_REGV_OFF_cit fundregv_it = regv_to_fundregv_off.find(regv);
        assert(fundregv_it != regv_to_fundregv_off.end());
        REGV_TO_TRIV_OFF_cit triv_it = regv_to_triv_off.find(fundregv_it->second.first);
        if (triv_it != regv_to_triv_off.end())
        {
          // store the mapping regv -> triv + offset
          regv_to_triv_off[regv] = std::make_pair(
            triv_it->second.first,
            fundregv_it->second.second);
        }
      }
    }


    TRIF_TO_REGF trif_to_regf;
    HASH_TO_TRIF hash_to_trif;
    REGF_TO_TRIF_CYCLE regf_to_trif;

    // FACES
    // Take only triangles which have a vertex in the fundamental domain.
    // According to DolHus96, see notes p.3299, these are periodic.
    // We shift them to be canonical triangles (note that we cannot filter only
    // the canonical ones, because they do not necessarily have a point in the
    // fundamental domain).
    for (Reg_Finite_faces_iterator regf_fit = reg.finite_faces_begin();
         regf_fit != reg.finite_faces_end();
         ++regf_fit)
    {
      Reg_Face_handle regf = regf_fit; // cast operator
      // get corresponding vertices and offsets in the periodic tri:
      Offset off[3];
      Del_Vertex_handle triv[3];
      bool one_vertex_in_domain = false;
      for (int i=0; i<3; i++)
      {
        REGV_TO_TRIV_OFF_cit triv_it = regv_to_triv_off.find(regf->vertex(i));
        // if there is a vertex image which is hidden in the fundamental domain,
        // then the triangle cannot be one of those we search
        if (triv_it == regv_to_triv_off.end())
        {
          one_vertex_in_domain = false;
          break;
        }
        triv[i] = triv_it->second.first;
        off[i] = triv_it->second.second;
        if (off[i].is_zero())
          one_vertex_in_domain = true;
      }

      if (not one_vertex_in_domain) continue;

      // set offsets: we might need to shift the face to have a canonical face
      int shift_x = std::min(off[0].x(), std::min(off[1].x(), off[2].x()));
      int shift_y = std::min(off[0].y(), std::min(off[1].y(), off[2].y()));
      assert(shift_x >= -1);
      assert(shift_y >= -1);
      Shift shift((shift_x<0 ? 1 : 0), (shift_y<0 ? 1 : 0));
      for (int i=0; i<3; i++) {
        off[i] += shift;
      }

      // create a face in tri
      // to identify equivalent triangles, we use a unique fingerprint
      size_t hash = del_unique_face_hash(triv, off, tri);
      if (hash_to_trif.find(hash) == hash_to_trif.end())
      {
        // keep the same order of vertices as in regf:
        Del_Face_handle trif = tri.tds().create_face(triv[0], triv[1], triv[2]);
        trif->set_offsets(tri.off_to_int(off[0]),
                          tri.off_to_int(off[1]),
                          tri.off_to_int(off[2]));
        assert(del_face_is_canonical(trif, tri));
        trif->info().dual_is_canonical = true;
        // set faces
        for (int i=0; i<3; i++) {
          triv[i]->set_face(trif);
        }
        // for setting neighbors, we need the map trif -> regf
        // order of vertices is the same in both faces
        trif_to_regf[trif] = regf;
        regf_to_trif[regf] = std::make_pair(trif, 0);
        hash_to_trif[hash] = trif;
      }
      else
      {
        // when we store the same trif for several regf, it might happen that
        // the non-original regf is turned, as compared to the original one
        // (which is stored as trif_to_regf)
        Del_Face_handle trif = hash_to_trif[hash];
        assert(trif_to_regf.find(trif) != trif_to_regf.end());
        Reg_Face_handle regf_orig = trif_to_regf[trif];
        assert(regv_to_fundregv_off.find(regf->vertex(0)) != regv_to_fundregv_off.end());
        Reg_Vertex_handle regv = regv_to_fundregv_off.find(regf->vertex(0))->second.first;

        int cycle;
        for (cycle=0; cycle<3; cycle++) {
          assert(regv_to_fundregv_off.find(regf_orig->vertex(cycle)) != regv_to_fundregv_off.end());
          Reg_Vertex_handle regv_orig = regv_to_fundregv_off.find(regf_orig->vertex(cycle))->second.first;
          if (regv == regv_orig)
            break; // keep this value of cycle
        }
        // thus, regf->vertex(i) == trif->vertex((i+cycle)%3)
        // because trif and regf_orig have same vertex numbering
        assert(cycle < 3); // make sure we found one
        regf_to_trif[regf] = std::make_pair(trif, cycle);
      }
    }

    // Set neighbors of the faces.
    // We loop over the same regular faces with one vertex in the fundamental
    // domain. Not all of their neighbors are periodic, but we will find all of
    // them if we loop over all such regular faces.
    for (REGF_TO_TRIF_CYCLE_it regf_it = regf_to_trif.begin();
         regf_it != regf_to_trif.end();
         ++regf_it)
    {
      Reg_Face_handle regf = regf_it->first;
      // vertices and neighbors are not necessarily in the same order in trif and regf:
      // Loop over the reg nbrs -- if it is periodic, we found also a tri nbr
      for (int i=0; i<3; i++)
      {
        Reg_Face_handle regnbr = regf->neighbor(i);
        if (regf_to_trif.find(regnbr) != regf_to_trif.end())
        {
          Del_Face_handle trif = regf_it->second.first;
          int j = (i + regf_it->second.second) % 3;

          Del_Face_handle trifnbr = regf_to_trif[regnbr].first;
          assert(trif->neighbor(j) == trifnbr or trif->neighbor(j) == Del_Face_handle(NULL));
          trif->set_neighbor(j, trifnbr);
        }
      }
    }

    # ifndef NDEBUG
    // Check that all faces have their neighbors
    for (Del_Face_iterator trif = tri.faces_begin();
         trif != tri.faces_end();
         ++trif)
    {
      for (int i=0; i<3; i++)
      {
        if (trif->neighbor(i) == Del_Face_handle(NULL))
        {
          Reg_Face_handle regf = trif_to_regf[trif];
          Reg_Face_handle regnbr = regf->neighbor(i);
          std::cerr << "ERROR: not all required neighbors are periodic in regular triangulation." << std::endl;
          std::cerr << "problematic face:"
                    << " " << &*trif << " (regf " << &*regf << ") nbr " << i << " (nbr regf " << &*regnbr << ")";
          for (int j=0; j<3; j++) {
            std::cerr << "  (" << regnbr->vertex(j)->point().point() << ") ";
            //std::cerr << "  (" << triv[j]->info().unique_id << " " << off[j].x() << " " << off[j].y() << ")";
          }
          std::cerr << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    # endif

    //del_debug_nbrs(tri, hash_to_trif);

  } // >>>
  void del_debug_nbrs(const Tri& tri, const HASH_TO_TRIF& hash_to_trif) const // <<<
  // This function helped to identify problems with CGAL::Handle_hash_function,
  // and with faces that have the same neighbor twice
  {
    # define PRINT(x) x
    del_save("del_debug.mytri", tri);

    // make sure that every face has neighbors
    PRINT(std::cout << "Listing nbrs:" << std::endl;)
    for (Del_Face_iterator trif = tri.faces_begin();
         trif != tri.faces_end();
         ++trif)
    {
      PRINT(std::cout << "trif=" << &*trif << std::endl;)
      for (int i=0; i<3; i++) {
        Del_Face_handle nbr = trif->neighbor(i);
        PRINT(std::cout << "   nbr=" << &*nbr << "  ";)
        PRINT(for (int j=0; j<3; j++) std::cout << " " << &*nbr->neighbor(j);)
        PRINT(std::cout << std::endl;)
      }
    }

    PRINT(std::cout << "Testing nbrs:" << std::endl;)
    for (Del_Face_iterator trif = tri.faces_begin();
         trif != tri.faces_end();
         ++trif)
    {
      PRINT(std::cout << "trif=" << &*trif << std::endl;)
      for (int i=0; i<3; i++) {
        assert(trif->neighbor(i) != Del_Face_handle(NULL));
        assert(trif->index(trif->neighbor(i)) == i);
        Del_Face_handle nbr = trif->neighbor(i);
        PRINT(std::cout << "   nbr=" << &*nbr << "  ";)
        PRINT(for (int j=0; j<3; j++) std::cout << " " << &*nbr->neighbor(j);)
        PRINT(std::cout << std::endl;)
        assert(nbr->has_neighbor(trif));
        # ifndef NDEBUG
        int n = nbr->index(trif);
        assert(0 <= n and n < 3);
        assert(trif->vertex(ccw(i)) == nbr->vertex(cw(n)));
        assert(trif->vertex(cw(i)) == nbr->vertex(ccw(n)));
        # endif
      }
    }

    for (typename std::map<size_t, Del_Face_handle>::const_iterator ut_it = hash_to_trif.cbegin();
         ut_it != hash_to_trif.cend();
         ++ut_it)
    {
      Del_Face_handle trif = ut_it->second;
      PRINT(std::cout << "trif=" << &*trif << std::endl;)
      for (int i=0; i<3; i++) {
        assert(trif->neighbor(i) != Del_Face_handle(NULL));
        assert(trif->index(trif->neighbor(i)) == i);
        Del_Face_handle nbr = trif->neighbor(i);
        PRINT(std::cout << "   nbr=" << &*nbr << "  ";)
        PRINT(for (int j=0; j<3; j++) std::cout << " " << &*nbr->neighbor(j);)
        PRINT(std::cout << std::endl;)
        # ifndef NDEBUG
        int n = nbr->index(trif);
        assert(0 <= n and n < 3);
        assert(trif->vertex(ccw(i)) == nbr->vertex(cw(n)));
        assert(trif->vertex(cw(i)) == nbr->vertex(ccw(n)));
        # endif
      }
    }

    # undef PRINT
  } // >>>
  void del_comply_internals(Tri& tri) const // <<<
  // Make tri comply with CGAL, concerning private and protected objects.
  // NOTE that we still do not pass the test tri.is_valid(), because we switch only to 3x3 if the 2-cycle test fails,
  // and is_valid() insists on the (stronger) test that there be no too_long_edges.
  {
    assert(tri.is_1_cover());

    // our domain is not square:
    tri._edge_length_threshold = FT(2*0.166)*CGAL::sqrt(_halfway_squared);

    // reset the internal variables
    tri._too_long_edges.clear();
    tri._too_long_edge_counter = 0;
    tri._virtual_vertices.clear();
    tri._virtual_vertices_reverse.clear();

    // build list of too long edges
    // as in Periodic_2_triangulation_2::load
    int i = 0;
    for (Del_Vertex_iterator vi = tri.vertices_begin();
         vi != tri.vertices_end(); ++vi)
    {
      tri._too_long_edges[vi] = std::list<Del_Vertex_handle>();
      ++i;
    }
    tri._too_long_edge_counter = tri.find_too_long_edges(tri._too_long_edges);
  } // >>>

public:
// simple functions concerning/using Delaunay structure
  Weight del_largest_weight(const Tri& tri) const // <<<
  {
    Weight max_weight(-1);
    for (Del_Unique_vertex_iterator dv_it = tri.unique_vertices_begin(),
                                    dv_end = tri.unique_vertices_end();
         dv_it != dv_end; ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      max_weight = std::max(max_weight, dv->info().weight);
    }
    assert(max_weight > 0);
    return max_weight;
  } // >>>
  Point weighted_circumcenter(Weight& weight, const Del_Face_handle& df, const Tri& tri, const Reg& reg = Reg()) const // <<<
  {
    // TODO make this nicer
    Weighted_point wp0(unwrapped(df->vertex(0)->point(), tri.get_offset(df, 0)), df->vertex(0)->info().weight);
    Weighted_point wp1(unwrapped(df->vertex(1)->point(), tri.get_offset(df, 1)), df->vertex(1)->info().weight);
    Weighted_point wp2(unwrapped(df->vertex(2)->point(), tri.get_offset(df, 2)), df->vertex(2)->info().weight);
    Weighted_point wp = reg.weighted_circumcenter(wp0, wp1, wp2);
    weight = squared_distance(wp.point(), df->vertex(0)->point()) - df->vertex(0)->info().weight;
    return wp.point();
  } // >>>

// functions concerning/using the Voronoi structure
  bool vor_vertex_is_covered_by_spheres(const Vor_vertex& vv, const Del_Face_handle& start, const Shift& start_shift, const Tri& tri, const Weight& max_radius_squared) const // <<<
  // In the vicinity of the Voronoi vertex search for spheres
  // that are closer than their radius.
  {
    // We keep a list of neighbors which are to be searched for close spheres.
    // (A PeriFace is a face + an offset shift which
    // puts the face into its right position relative to the
    // original vertex vv)
    typedef std::pair<Del_Face_handle, Shift> PeriFace;
    typedef typename std::unordered_set<PeriFace, HOhash> PeriFaces;
    typedef typename PeriFaces::iterator PeriFaces_iterator;
    PeriFaces dfaces, dfaces_tocheck, dfaces_checked;

    // Initialise dfaces:
    // We did not use vv.offset when finding start!
    dfaces.insert(std::make_pair(start, start_shift));

    // repeat to add neighbors of faces until we are done
    while (not dfaces.empty())
    {
      dfaces_tocheck.clear();
      for (PeriFaces_iterator df_it = dfaces.begin();
           df_it != dfaces.end();
           ++df_it)
      {
        const Del_Face_handle df = df_it->first;;
        const Shift shift = df_it->second;
        bool must_add_nbrs = false;

        for (int i=0; i<3; i++) {
          Point p = unwrapped(df->vertex(i)->point(), tri.get_offset(df, i) + shift);
          // check the distance between vv and Delaunay vertices (spheres)
          if (squared_distance(vv.point, p) < df->vertex(i)->info().weight)
            return true;
          // if we are closer than max_radius, add all nbrs
          if (squared_distance(vv.point, p) < max_radius_squared)
            must_add_nbrs = true;
        }

        if (must_add_nbrs) {
          for (int i=0; i<3; i++)
            dfaces_tocheck.insert(std::make_pair(
                df->neighbor(i),
                tri.combine_offsets(shift, -tri.get_neighbor_offset(df, i)))); // XXX
        }
        dfaces_checked.insert(*df_it);
      } // end loop over dfaces

      // prepare the new dfaces from the dfaces_tocheck
      dfaces.clear();
      for (PeriFaces_iterator df_it = dfaces_tocheck.begin();
           df_it != dfaces_tocheck.end();
           ++df_it)
      {
        if (dfaces_checked.find(*df_it) == dfaces_checked.end())
          dfaces.insert(*df_it);
      }
    } // end while

    // we did not find any overlapping sphere
    return false;
  } // >>>
  bool vor_edge_is_split_by_sphere(const Point& p0, const Point& p1, const Point& center, const Weight& radius_sq) const // <<<
  // TODO: make this a predicate of a Geometry_traits
  // -- just in case we want to use exact constructions at some point ...
  // XXX We may assume here that the two boundary points are not covered.
  {
    Compute_scalar_product csp;
    Construct_perpendicular_vector cpv;
    Vector parall(p0, p1), // p1 - p0
           vec0(center, p0), vec1(center, p1);
    Vector perpen = cpv(parall, CGAL::COUNTERCLOCKWISE);
    FT proj0 = csp(parall, vec0),
       proj1 = csp(parall, vec1),
       ortho = csp(perpen, vec0); // signed distance center--straight
    FT tmp = radius_sq*csp(parall, parall) - ortho*ortho;

    // if one of p0, p1 is covered: segment is covered as well
    if (proj0*proj0 <= tmp or proj1*proj1 <= tmp) return true;

    // p0, p1 are not within the sphere. If they are on two different sides of
    // it, the sphere can cover the segment.
    if (proj0*proj1 < 0 and tmp > 0) return true;
    return false;
  } // >>>
  bool vor_face_contains_point(const Point& p, const Del_Vertex_handle& dv, const Offset& dv_off, const Tri& tri) const // <<<
  // dv_off is the offset of dv, plus applied all shifts we need to put it locally around p
  // Moreover, it avoids repeated evaluations of tri.get_offset(dv). (See discussion around vor_locate, notes p.3302)
  // NOTE that this function is not stable against degenerate cases. Use
  // vor_vertex_sector_of_point if you really must rely on the result.
  {
    // circulate around the Delaunay vertex to get all Voronoi edges
    // then ask whether p is on the right side of all of them
    const Del_My_Face_circulator df0_started(dv);
    Del_My_Face_circulator df0_it(df0_started), df1_it(df0_started);
    df1_it++; // df1_it is one face ahead of df0_it
    do {
      const Vor_vertex& vv0 = df0_it.base()->info();
      const Vor_vertex& vv1 = df1_it.base()->info();
      // vv.point + domain*vv.offset   is location w.r.t. the dual DelFace of vv
      // we must get a shift to put it locally around dv
      Shift sh0 = dv_off - tri.get_offset(df0_it.base(), df0_it.index());
      Shift sh1 = dv_off - tri.get_offset(df1_it.base(), df1_it.index());

      CGAL::Orientation ori = tri.orientation(
        unwrapped(vv0.point, vv0.offset + sh0),
        unwrapped(vv1.point, vv1.offset + sh1),
        p);
      if (ori == CGAL::RIGHT_TURN)
        return false;

      // increment
      ++df1_it, ++df0_it;
    } while (df0_it != df0_started);
    return true;
  } // >>>
  bool vor_face_is_convex(const Del_Vertex_handle& dv, const Tri& tri) const // <<<
  {
    // Circulate with three circulators over vertices of the VorFace(dv)
    const Del_My_Face_circulator df_started(dv);
    Del_My_Face_circulator df_upp_it(df_started), df_mid_it(df_started), df_low_it(df_started);
    df_upp_it++; df_upp_it++;
    df_mid_it++;

    do {
      const Vor_vertex& vv_low = df_low_it.base()->info(),
                        vv_mid = df_mid_it.base()->info(),
                        vv_upp = df_upp_it.base()->info();
      Shift shift_mid = tri.get_offset(df_low_it.base(), df_low_it.index()) - tri.get_offset(df_mid_it.base(), df_mid_it.index()),
            shift_upp = tri.get_offset(df_low_it.base(), df_low_it.index()) - tri.get_offset(df_upp_it.base(), df_upp_it.index());
      Point vp_low = unwrapped(vv_low.point, vv_low.offset),
            vp_mid = unwrapped(vv_mid.point, vv_mid.offset + shift_mid),
            vp_upp = unwrapped(vv_upp.point, vv_upp.offset + shift_upp);
      if (tri.orientation(vp_low, vp_mid, vp_upp) == CGAL::RIGHT_TURN)
        return false;
      ++df_upp_it, ++df_mid_it, ++df_low_it;
    } while (df_low_it != df_started);
    return true;
  } // >>>
  bool vor_faces_are_convex(const Tri& tri) const // <<<
  {
    for (Del_Unique_vertex_iterator dv_it = tri.unique_vertices_begin(),
                                    dv_end = tri.unique_vertices_end();
         dv_it != dv_end; ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      if (not vor_face_is_convex(dv, tri))
        return false;
    }
    return true;
  } // >>>

// functions which must have a clear answer, independent on precision (predicates?)
  int vor_vertex_sector_of_point(const Del_Face_handle& df, const Shift& df_shift, const Point& point, const Tri& tri) const // <<<
  // For the given VorVertex(df+df_shift) returns the index of the
  // DelVertex(VorFace) in which is the point. This question must have an
  // answer, even if the point is on a VorEdge, or on the VorVertex itself.
  //
  // NOTE: that the answer is not always the same in degenerate cases.
  //
  // Algorithm:
  //   a Loop over edges and store the orientation decisions.
  //   b Return the VorFacet corresponding to the found sector.
  //   c If no sector found (all three edges give the same answer)
  //     the point must be enclosed in a precision-size triangle where
  //     all three edges meet.
   {
    const Vor_vertex& vv = df->info();
    Point pv = unwrapped(vv.point, vv.offset + df_shift);

    // a
    CGAL::Orientation orients[3];
    for (int i=0; i<3; i++) {
      Del_Face_handle df_nbr = df->neighbor(i);
      const Vor_vertex& vv_nbr = df_nbr->info();
      Shift shift = tri.combine_offsets(df_shift, -tri.get_neighbor_offset(df, i));
      Point pv_nbr = unwrapped(vv_nbr.point, vv_nbr.offset + shift);
      orients[i] = tri.orientation(pv, pv_nbr, point);
    }

    // b
    int take_dv_i = -1;
    for (int i=0; i<3; i++) {
      if (orients[i] != CGAL::RIGHT_TURN and orients[ccw(i)] != CGAL::LEFT_TURN)
      {
        // In a degenerate case, several sectos may claim the point
        // Just take the first...
        take_dv_i = cw(i); // see notes p.3307 for a picture
        break;
      }
    }

    // c
    if (take_dv_i == -1) {
      std::cerr << "INFO: Three Voronoi edges do not meet in a point but in a small triangle."
                << " Dist_sq=" << squared_distance(point, pv)
                << std::endl;
      assert(squared_distance(point, pv) < PREC);
      take_dv_i = 0; // just take any facet
    }

    return take_dv_i;
  } // >>>
  bool vor_face_sector_contains_point(const Del_Face_handle& df, const Shift& df_shift, const int dv_i, const Point& point, const Clusters& clusters, DF_TO_CLS& df_to_cls, const Tri& tri) const // <<<
  // Precondition: The point is in the VorFace(dv), or at least in the sector as seen from VorVertex(df).
  //
  // Return value says whether the point is in that part of the VorFace which
  // belongs to VorVertex (as seen from dv). // We do not check overlap with sphere here...
  // This question has a precise answer, because
  //   either the whole edge belongs to the same cluster
  //   or the sphere cuts the edge an then precision does not play a role
  //
  // Algorithm (see notes p.3308 for a picture)
  // VorVertices: vv_low -> me -> vv_upp
  // 1. If the point is not in the sector (vv_low, vv_upp),
  //    it cannot be in the part of vv
  // 2. The point is either in the sector (vv_low,vv) or in (vv,vv_upp)
  //    Determine in which.
  // 3. Test the over vv for same cluster
  // 4. If different cluster, construct perpendicular drop
  //    and decide by geometry.
  //
  {
    // me
    const Vor_vertex& vv = df->info();
    Point vp = unwrapped(vv.point, vv.offset + df_shift);

    // sphere
    Del_Vertex_handle dv = df->vertex(dv_i);
    Offset dv_off = tri.get_offset(df, dv_i) + df_shift;
    Point dp = unwrapped(dv->point(), dv_off); // can be in the VorFace(dv) or not!

    // edge-neighbors
    Del_My_Face_circulator df_low_it(dv, df, dv_i), df_upp_it(dv, df, dv_i);
    --df_low_it, ++df_upp_it;
    const Vor_vertex& vv_low = df_low_it.base()->info();
    const Vor_vertex& vv_upp = df_upp_it.base()->info();
    Shift vv_low_shift = dv_off - tri.get_offset(df_low_it.base(), df_low_it.index()),
          vv_upp_shift = dv_off - tri.get_offset(df_upp_it.base(), df_upp_it.index());
    Point vp_low = unwrapped(vv_low.point, vv_low.offset + vv_low_shift),
          vp_upp = unwrapped(vv_upp.point, vv_upp.offset + vv_upp_shift);
    assert(tri.orientation(vp_low, vp, vp_upp) == CGAL::LEFT_TURN);

    // we need some differentiation to allow dp not in the face
    Del_Face_handle df_other;
    Point vp_other;
    if (tri.orientation(vp_low, vp, dp) == CGAL::LEFT_TURN)
    { // the triangle (dp, vp_low, vp) is well-oriented (positive face sign)

      if (tri.orientation(vp, vp_upp, dp) == CGAL::LEFT_TURN)
      { // the triangle (dp, vp, vp_upp) is well-oriented (positive face sign)

        // test the wider sector (vp_low, vp_upp):
        if (tri.orientation(dp, vp_low, point) == CGAL::RIGHT_TURN or
            tri.orientation(dp, vp_upp, point) == CGAL::LEFT_TURN)
          return false;

        // choose vp_low or vp_upp
        // work with df_upp from now on
        if (tri.orientation(dp, vp, point) == tri.orientation(dp, vp, vp_low)) {
          df_other = df_low_it.base();
          vp_other = vp_low;
        } else {
          df_other = df_upp_it.base();
          vp_other = vp_upp;
        }
      } else
      { // the triangle (dp, vp, vp_upp) is to be ignored

        // test the wider sector (vp_low, vp)
        if (tri.orientation(dp, vp_low, point) == CGAL::RIGHT_TURN or
            tri.orientation(dp, vp, point) == CGAL::LEFT_TURN)
          return false;

        // use vp_low for the perpendicular drop
        df_other = df_low_it.base();
        vp_other = vp_low;
      }
    } else
    { // the triangle (dp, vp_low, vp) is to be ignored

      if (tri.orientation(vp, vp_upp, dp) == CGAL::LEFT_TURN)
      { // the triangle (dp, vp, vp_upp) is well-oriented (positive face sign)

        // test the wider sector (vp_upp, vp)
        if (tri.orientation(dp, vp, point) == CGAL::RIGHT_TURN or
            tri.orientation(dp, vp_low, point) == CGAL::LEFT_TURN)
          return false;

        // use vp_low for the perpendicular drop
        df_other = df_upp_it.base();
        vp_other = vp_upp;
      } else
      { // the triangle (dp, vp, vp_upp) is to be ignored
        return false;
      }
    }

    // test whether the other end of the edge is in a cluster
    if (df_other->info().is_covered)
      return true;

    // no problem if the whole edge is in the cluster
    // XXX this is not sufficient to exclude a split edge!
    if (clusters.find(df_to_cls[df]) == clusters.find(df_to_cls[df_other]))
      return true;

    // the basins of "responsibility" of the two Voronoi vertices
    // are separated by the perpendicular drop on the edge
    Construct_perpendicular_vector cpv;
    Compute_scalar_product csp;
    Vector edge(vp, vp_other), // vp_other - vp
           perp(cpv(edge, CGAL::COUNTERCLOCKWISE));
    FT tmp = csp(perp, Vector(dp, vp)) / csp(perp, perp);
    Point E(dp.x() + tmp*perp.x(),
            dp.y() + tmp*perp.y());
    return (tri.orientation(dp, E, point) == tri.orientation(dp, E, vp));
  } // >>>
  Del_My_Face_circulator vor_face_sector_containing_point(const Del_Vertex_handle& dv, const Offset& dv_off, const Point& point, const Tri& tri, const FT& prec=-1.0) const // <<<
  // Returns the Voronoi vertex (+index of dv) in whose basin of responsability is the point.
  // CutCells with negative face sign must be treated specially..
  // The optional precision parameter allows to favour non-covered vertices.
  //
  // Precondition: We assume that the point is in the Voronoi face!
  //
  // Algorithm: see notes p.3309 for a picture
  // 1. Identify in which triangle (dv, vv_low, vv_upp) the point is.
  //    This steps throws away the badly oriented triangles
  //    There must be one such triangle (or, there is a negative facet sign,
  //    and the point is identical to the boundary VorVertices)
  // 2. It remains to choose between two points.
  //    Construct the perpendicular drop and decide (just return one of them if
  //    they are on the same side)
  //    In this last step, we might want to avoid to return a "covered" VorVertex.
  //    (Need a precision parameter for this)
  {

    const Point dp = unwrapped(dv->point(), dv_off);
    Del_My_Face_circulator df_it_mid, df_it_upp;
    Point vp_mid, vp_upp;
    { // Circulate around dv and identify the AVV triangle in which is the point.
      # ifndef NDEBUG
      bool have_avv = false;
      # endif
      const Del_My_Face_circulator df_started(dv);
      Del_My_Face_circulator df_it_low(df_started);
      df_it_mid = df_started; ++df_it_mid;
      df_it_upp = df_started; ++df_it_upp, ++df_it_upp;
      Point
      vp_low = unwrapped(df_it_low.base()->info().point,
                         df_it_low.base()->info().offset + dv_off - tri.get_offset(df_it_low.base(), df_it_low.index())); // place df (and vv) with respect to dv
      vp_mid = unwrapped(df_it_mid.base()->info().point,
                         df_it_mid.base()->info().offset + dv_off - tri.get_offset(df_it_mid.base(), df_it_mid.index())); // place df (and vv) with respect to dv
      CGAL::Orientation ori_low = tri.orientation(dp, vp_low, point),
                        ori_mid = tri.orientation(dp, vp_mid, point);
      bool prev_was_inverted = (tri.orientation(vp_low, vp_mid, dp) == CGAL::RIGHT_TURN);

      do {
        vp_upp = unwrapped(df_it_upp.base()->info().point,
                           df_it_upp.base()->info().offset + dv_off - tri.get_offset(df_it_upp.base(), df_it_upp.index())); // place df (and vv) with respect to dv
        CGAL::Orientation ori_upp = tri.orientation(dp, vp_upp, point);

        // It may happen that the Delaunay vertex (dp) is not inside the Voronoi facet.
        // We consider the point to be in one of the triangles *only* if the triangle is not "inverted".
        bool is_inverted = (tri.orientation(vp_mid, vp_upp, dp) == CGAL::RIGHT_TURN);
        if (is_inverted)
        { // maybe the last triangle was OK
          if ((not prev_was_inverted) and ori_low != CGAL::RIGHT_TURN) {
            // the last triangle was degenerate, and the point was close to the boundary
            // turn back the clock:
            df_it_upp = df_it_mid; vp_upp = vp_mid;
            df_it_mid = df_it_low; vp_mid = vp_low;
            # ifndef NDEBUG
            have_avv = true;
            # endif
            break;
          }
        } else if (ori_upp != CGAL::LEFT_TURN and (ori_mid != CGAL::RIGHT_TURN or prev_was_inverted))
        {
          # ifndef NDEBUG
          have_avv = true; // we found the AVV triangle
          # endif
          break;
        }

        // increment
        prev_was_inverted = is_inverted;
        ori_low = ori_mid, vp_low = vp_mid;
        ori_mid = ori_upp, vp_mid = vp_upp;
        ++df_it_low, ++df_it_mid, ++df_it_upp;
      } while (df_it_low != df_started);
      // Must succeed
      assert(have_avv);
    }

    // Identify which of the two vv has the correct cluster:
    // (Must also succeed)
    // the basins of "responsibility" of the two Voronoi vertices
    // are separated by the perpendicular drop on the edge
    Construct_perpendicular_vector cpv;
    Compute_scalar_product csp;
    Vector edge(vp_mid, vp_upp), // vp_upp - vp_mid
           perp(cpv(edge, CGAL::CLOCKWISE));
    FT tmp = csp(perp, Vector(dp, vp_mid)) / csp(perp, perp);
    Point E(dp.x() + tmp*perp.x(),
            dp.y() + tmp*perp.y());
    CGAL::Orientation
    ori_mid = tri.orientation(dp, E, vp_mid),
    ori_upp = tri.orientation(dp, E, vp_upp),
    ori_p   = tri.orientation(dp, E, point);
    // we want either left or right
    if (ori_mid == CGAL::COLLINEAR) ori_mid = CGAL::LEFT_TURN;
    if (ori_upp == CGAL::COLLINEAR) ori_upp = CGAL::LEFT_TURN;
    if (ori_p   == CGAL::COLLINEAR) ori_p   = CGAL::LEFT_TURN;

    Del_My_Face_circulator* result;
    if (ori_mid != ori_upp)
    { // The perpendicular drop E is between the two Voronoi vertices

      if (ori_p == ori_mid) result = &df_it_mid;
      else result = &df_it_upp;

      // avoid to return a covered vertex
      // make an exception with precision prec
      if (prec > 0) {
        if (result->base()->info().is_covered and CGAL::abs(csp(edge, Vector(point, E))) < prec)
          // take the other one:
          result = (result == &df_it_mid ? &df_it_upp : &df_it_mid);
      }
    } else {
      result = (df_it_mid.base()->info().is_covered ? &df_it_upp : &df_it_mid);
    }

    return *result;
  } // >>>
  Del_My_Face_circulator del_vertex_sector_containing_point(const Del_Vertex_handle& dv, const Offset& dv_off, const Point& point, const Tri& tri, __attribute__((unused)) const FT& prec=-1.0) const // <<<
  // For the given DelVertex (+offset) returns that sector of space
  // (DelFace = VorVertex) which contains the point.
  //
  // This function can replace vor_face_sector_containing_point
  // Precondition: We assume that the point is in the Voronoi face!
  //
  // see notes p.3309 for a picture
  // Algorithm: as vor_vertex_sector_of_point
  //   a Loop over edges and store the orientation decisions.
  //   b Return the DelFacet corresponding to the found sector.
  //     (exactly one has 2 hits)
  //   c If no sector found, the case is degenerate (point close to dv)
  {
    const Point dp = unwrapped(dv->point(), dv_off);

    // knows which vertex is dv while circulating around it
    const Del_My_Face_circulator df_started(dv);
    Del_My_Face_circulator df_it(df_started);

    Del_Face_handle df = df_it.base();
    int i = cw(df_it.index()); // test the "forward edge" of the DelFace
    Shift sh = dv_off - tri.get_offset(df, df_it.index()); // the shift which must be applied to df in order to put it locally around the correct copy of dv
    Offset off = tri.get_offset(df, i) + sh; // the offset of the vertex we test, to put it locally around dv
    Point nbr = unwrapped(df->vertex(i)->point(), off);
    CGAL::Orientation prev_ori = tri.orientation(dp, nbr, point);
    # ifndef NDEBUG
    CGAL::Orientation first_ori = prev_ori;
    # endif
    ++df_it;
    while (df_it != df_started)
    {
      df = df_it.base();
      i = cw(df_it.index()); // test the "forward edge" of the DelFace
      sh = dv_off - tri.get_offset(df, df_it.index());
      off = tri.get_offset(df, i) + sh;
      nbr = unwrapped(df->vertex(i)->point(), off);
      CGAL::Orientation ori = tri.orientation(dp, nbr, point);
      if (prev_ori == CGAL::LEFT_TURN and ori == CGAL::RIGHT_TURN)
        return df_it;
      prev_ori = ori;
      ++df_it;
    }
    // Either the point is in the starting df, or it is degenerate: in both cases the starting df can be returned
    # ifndef NDEBUG
    // test the starting df as well:
    if (prev_ori == CGAL::LEFT_TURN and first_ori == CGAL::RIGHT_TURN)
      return df_it; // in starting df
    assert(squared_distance(dp, point) < prec); // degenerate
    # endif
    return df_it;
  } // >>>

// functions for periodicity (independent of all weights and triangulation types)
  Point wrapped(Offset& off, const Point& p) const // <<<
  // similar to Periodic_triangulation: move_in_domain
  //   p = result + domain*off
  {
    FT x = p.x();
    FT y = p.y();

    // zero the offset:
    off -= off;

    while (x <  _domain.xmin()) { x += _domain.xmax() - _domain.xmin(); off.x()--; }
    while (x >= _domain.xmax()) { x -= _domain.xmax() - _domain.xmin(); off.x()++; }

    while (y <  _domain.ymin()) { y += _domain.ymax() - _domain.ymin(); off.y()--; }
    while (y >= _domain.ymax()) { y -= _domain.ymax() - _domain.ymin(); off.y()++; }

    return Point(x, y);
  } // >>>
  Point wrapped(const Point& p, const FT& eps=-1.0) const // <<<
  // same as Periodic_triangulation: move_in_domain
  {
    FT x = p.x();
    FT y = p.y();

    while (x <  _domain.xmin()) x += _domain.xmax() - _domain.xmin();
    while (x >= _domain.xmax()) x -= _domain.xmax() - _domain.xmin();
    // it may happen that x still lies outside ...
    if (eps > 0 and x <= _domain.xmin()) x = _domain.xmin() + eps;

    while (y <  _domain.ymin()) y += _domain.ymax() - _domain.ymin();
    while (y >= _domain.ymax()) y -= _domain.ymax() - _domain.ymin();
    if (eps > 0 and y <= _domain.ymin()) y = _domain.ymin() + eps;

    return Point(x, y);
  } // >>>
  Point unwrapped(const Point& p, const Offset& o) const // <<<
  {
    return Point(p.x() + (_domain.xmax() - _domain.xmin()) * o.x(),
                 p.y() + (_domain.ymax() - _domain.ymin()) * o.y());
  } // >>>
  bool del_face_is_canonical(const Del_Face_handle& df, const Tri& tri) const // <<<
  // this code is copied from Periodic_2_triangulation_triangle_iterator_2 at
  // Periodic_2_triangulation_iterators_2.h
  {
    // fetch all offsets
    //Offset off0, off1, off2;
    //get_edge_offsets(off0, off1, off2);
    Offset face_off0 = tri.int_to_off(df->offset(0));
    Offset face_off1 = tri.int_to_off(df->offset(1));
    Offset face_off2 = tri.int_to_off(df->offset(2));
    Offset diff_off((   face_off0.x() == 1
                     && face_off1.x() == 1
                     && face_off2.x() == 1) ? -1 : 0,
                    (   face_off0.y() == 1
                     && face_off1.y() == 1
                     && face_off2.y() == 1) ? -1 : 0);
    Offset off0 = tri.combine_offsets(tri.get_offset(df, 0), diff_off);
    Offset off1 = tri.combine_offsets(tri.get_offset(df, 1), diff_off);
    Offset off2 = tri.combine_offsets(tri.get_offset(df, 2), diff_off);

    if (not tri.is_1_cover())
      {
        // If there is one offset with entries larger than 1 then we are
        // talking about a vertex that is too far away from the original
        // domain to belong to a canonical triangle.
        if (off0.x() > 1) return false;
        if (off0.y() > 1) return false;
        if (off1.x() > 1) return false;
        if (off1.y() > 1) return false;
        if (off2.x() > 1) return false;
        if (off2.y() > 1) return false;
      }

    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x() & off2.x();
    int offy = off0.y() & off1.y() & off2.y();

    return (offx == 0 && offy == 0);
  } // >>>
  size_t del_unique_face_hash(const Del_Face_handle& df, const Tri& tri) const // <<<
  {
    Del_Vertex_handle vertices[3];
    Offset offsets[3];
    for (int i=0; i<3; i++) {
      vertices[i] = tri.get_original_vertex(df->vertex(i));
      offsets[i] = tri.get_offset(df, i);
    }
    return del_unique_face_hash(vertices, offsets, tri);
  } // >>>
  size_t del_unique_face_hash(const Del_Vertex_handle* vertices, const Offset* offs, __attribute__((unused)) const Tri& tri) const // <<<
  // XXX make sure that the vertices are fundamental ones!
  {
    // we take the hash values of all three vertices
    // and two relative offsets
    //
    // It may happen that two periodic images of triangles use the vertices in
    // different orders, so it is important that the order of the vertices is
    // not seen in the hash value. Also, we need to assert that the two
    // relative offsets are the same for both triangles.
    // We first sort the vertices by their individual hash value, and then take
    // the offsets relative to the first of these.
    // NOTE: we can get *identical* vertices, so we have to do a hierarchical sort
    size_t i=0, j=1, k=2;

    // sort by
    //   1. vertex id
    //   2. if identical ids, by offset.x()
    //   3. if these also identical, by offset.y()
    # define CMP(n,m) \
        vertices[n]->info().unique_id < vertices[m]->info().unique_id or        \
        (vertices[n]->info().unique_id == vertices[m]->info().unique_id and     \
         (offs[n].x() < offs[m].x() or                                          \
          (offs[n].x() == offs[m].x() and offs[n].y() < offs[m].y())))
    if (CMP(k,j)) std::swap(k, j);
    if (CMP(j,i)) std::swap(j, i);
    if (CMP(k,j)) std::swap(k, j);
    # undef CMP

    # if HASH_WITH_IDS
    // This is a non-probabilistic version, but limited in the number of vertices:
    // Put the hash values into different blocks of the 64bits of size_t:
    size_t n=16, m=4, O=7;
    #   ifndef NDEBUG
    const size_t N=(1<<n);
    const size_t M=(1<<m); // allows shifts from -7 to 8 (XXX we do not assert this);
    assert(3*n + 4*m <= 64);
    assert(tri.number_of_vertices() < N);
    for (int l=0; l<3; l++) assert(vertices[l]->info().unique_id >= 0 and vertices[l]->info().unique_id < N);
    assert(offs[j].x()-offs[i].x()+O >= 0 and offs[j].x()-offs[i].x()+O < M);
    assert(offs[j].y()-offs[i].y()+O >= 0 and offs[j].y()-offs[i].y()+O < M);
    assert(offs[k].x()-offs[i].x()+O >= 0 and offs[k].x()-offs[i].x()+O < M);
    assert(offs[k].y()-offs[i].y()+O >= 0 and offs[k].y()-offs[i].y()+O < M);
    #   endif // NDEBUG
    size_t result = (vertices[i]->info().unique_id
                  + (vertices[j]->info().unique_id << n)
                  + (vertices[k]->info().unique_id << (2*n))
                  + (size_t(offs[j].x() - offs[i].x() + O) << (3*n))
                  + (size_t(offs[j].y() - offs[i].y() + O) << (3*n+m))
                  + (size_t(offs[k].x() - offs[i].x() + O) << (3*n+2*m))
                  + (size_t(offs[k].y() - offs[i].y() + O) << (3*n+3*m)));
    //std::cout << std::endl;
    //std::cout << std::bitset<64>(vertices[i]->info().unique_id) << std::endl;
    //std::cout << std::bitset<64>(vertices[j]->info().unique_id << n) << std::endl;
    //std::cout << std::bitset<64>(vertices[k]->info().unique_id << (2*n)) << std::endl;
    //std::cout << std::bitset<64>(size_t(offs[j].x() - offs[i].x() + O) << (3*n)) << std::endl;
    //std::cout << std::bitset<64>(size_t(offs[j].y() - offs[i].y() + O) << (3*n+m)) << std::endl;
    //std::cout << std::bitset<64>(size_t(offs[k].x() - offs[i].x() + O) << (3*n+2*m)) << std::endl;
    //std::cout << std::bitset<64>(size_t(offs[k].y() - offs[i].y() + O) << (3*n+3*m)) << std::endl;
    //std::cout << "----------------------------------------------------------------" << std::endl;
    //std::cout << std::bitset<64>(result) << std::endl;
    return result;
    # else
    // combine the hash values:
    // this is the implementation of boost::combine_hash
    std::hash<int> int_hash;
    size_t result = 0;
    result ^= vertices[i]->info().unique_id + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= vertices[j]->info().unique_id + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= vertices[k]->info().unique_id + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(offs[j].x() - offs[i].x()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(offs[j].y() - offs[i].y()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(offs[k].x() - offs[i].x()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(offs[k].y() - offs[i].y()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    return result;
    # endif
  } // >>>
  bool contains_2_cycle(const Tri& tri) const // <<<
  // expensive test whether tri contains a 2-cycle
  // this is similar to Periodic_2_triangulation_2::is_triangulation_in_1_sheet
  // but without the cover_1 shortcut test
  // Note: we cannot rely on internal checks of tri
  //       as long as its tri._edge_length_threshold relies on a square domain
  //       and we cannot escape from its cover_1 tests
  {
    // see Caroli+Teillaud(2009), or notes p.3293
    if (all_edges_shorter_than_halfway(tri))
      return false;

    for (Del_Vertex_iterator dv_it = tri.vertices_begin();
         dv_it != tri.vertices_end();
         ++dv_it)
    {
      Del_Vertex_handle dv = dv_it; // cast operator

      std::set<Del_Vertex_handle> nbrs;
      size_t degree = 0;

      // circulate over nbr vertices (via incident faces)
      Del_My_Face_circulator df_started(dv);
      Del_My_Face_circulator df_it(df_started);
      do {
        nbrs.insert(df_it.base()->vertex(ccw(df_it.index())));
        degree++;
      } while (++df_it != df_started);
      if (degree != nbrs.size())
        return false;
    }
    return true;
  } // >>>
  bool all_edges_shorter_than_halfway(const Tri& tri) const // <<<
  // This reimplements Periodic_2_triangulation_2::is_extensible_triangulation_in_1_sheet_h1
  // without the cover_1 check
  {
    typedef typename Tri::Geom_traits::FT                 FT;
    typedef typename Tri::Segment                         Segment;
    typedef typename Tri::Periodic_segment                Periodic_Segment;
    typedef typename Tri::Periodic_segment_iterator       Periodic_segment_iterator;

    assert(_halfway_squared > 0);
    FT halfway_sq = _halfway_squared * tri.number_of_sheets()[0]*tri.number_of_sheets()[0]; // hx==hy
    Segment s;
    for (Periodic_segment_iterator psit = tri.periodic_segments_begin(Tri::UNIQUE);
         psit != tri.periodic_segments_end(Tri::UNIQUE);
         ++psit)
    {
      const Periodic_Segment& seg = *psit;
      //s = tri.construct_segment(seg); // XXX why is this protected ?
      s = tri.geom_traits().construct_segment_2_object()(
           seg[0].first, seg[1].first, seg[0].second, seg[1].second);
      if (s.squared_length() > halfway_sq)
        return false;
    }
    return true;
  } // >>>

// input, output
  void del_save(const char* filename, const Tri& tri, const bool verbose=true) const // <<<
  {
    if (tri.empty()) return;
    std::ofstream os(filename, std::ios::out);

    // write domain and cover sheets
    os << tri.domain() << std::endl
       << tri.number_of_sheets()[0] << " " << tri.number_of_sheets()[1] << std::endl;

    // write the vertices
    os << tri.number_of_vertices() << std::endl;
    for (Del_Vertex_iterator dv_it = tri.vertices_begin();
         dv_it != tri.vertices_end();
         ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      const Point p = dv_it->point();
      os << unwrapped(dv->point(), tri.get_offset(dv)) // NOT !!! dv_it->offset();
         << " " << dv->info().weight
         //<< " " << &*dv
         << " " << dv->info().unique_id
         << std::endl;
    }
    os << std::endl;

    // write faces
    os << tri.number_of_faces() << std::endl;
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it;
      for (int i=0; i<3; i++) {
        os << unwrapped(df->vertex(i)->point(), tri.get_offset(df, i)) << " ";
      }
      os << &*df
         //<< " " << del_unique_face_hash(df, tri)
         << " " << &*df->neighbor(0) << "_" << &*df->neighbor(1) << "_" << &*df->neighbor(2)
         //<< " " << &*df->vertex(0) << "_" << &*df->vertex(1) << "_" << &*df->vertex(2)
         << " " << df->vertex(0)->info().unique_id << "_" << df->vertex(1)->info().unique_id << "_" << df->vertex(2)->info().unique_id
         << std::endl;
    }
    os << std::endl;

    // we do not write neighbors

    os.close();
    if (verbose) {
      std::cerr << "INFO: Periodic weighted Delaunay triangulation written to file \"" << filename << "\"" << std::endl;
    }
  } // >>>
  void vor_save(const char* filename, const Tri& tri, const bool verbose=true) const // <<<
  {
    if (tri.empty()) return;
    std::ofstream os(filename, std::ios::out);

    // write domain and cover sheets
    os << tri.domain() << std::endl
       << tri.number_of_sheets()[0] << " " << tri.number_of_sheets()[1] << std::endl;

    // write the Voronoi vertices
    os << tri.number_of_faces() << std::endl;
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it;
      const Vor_vertex& vv = df->info();
      os << unwrapped(vv.point, vv.offset)
         << " " << vv.weight
         << " " << &*df
         << " " << vv.is_covered
         //<< " " << vv.cluster_id
         << std::endl;
    }
    os << std::endl;

    // write the Voronoi faces
    os << tri.number_of_vertices() << std::endl;
    int cnt = 0;
    for (Del_Vertex_iterator dv_it = tri.vertices_begin();
         dv_it != tri.vertices_end();
         ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      const Vor_face& vf = dv->info();

      const Del_My_Face_circulator df_started(dv);
      Del_My_Face_circulator df_it(df_started);
      do {
        Del_Face_handle df = df_it; // cast operator
        const Vor_vertex& vv = df->info();
        Shift shift = tri.get_offset(dv) - tri.get_offset(df, df_it.index());
        os << unwrapped(vv.point, vv.offset + shift) << " ";
      } while (++df_it != df_started);
      os << vf.weight
         << " " << &*dv
         << std::endl;
      cnt++;
    }
    os << std::endl;
    assert(cnt == tri.number_of_vertices());

    // we do not write neighbors.

    os.close();
    if (verbose) {
      std::cerr << "INFO: Periodic weighted Voronoi 'triangulation' written to file \"" << filename << "\"" << std::endl;
    }
  } // >>>
  void reg_save(const char* filename, const Reg& reg, const bool verbose=true) const // <<<
  {
    typedef typename Reg::Finite_vertices_iterator        Reg_Finite_vertices_iterator;
    typedef typename Reg::Finite_faces_iterator           Reg_Finite_faces_iterator;

    std::ofstream os(filename, std::ios::out);

    // write domain and cover sheets
    os << _domain << std::endl
       << 1 << " " << 1 << std::endl;

    // write the vertices
    os << reg.number_of_vertices() << std::endl;
    for (Reg_Finite_vertices_iterator dv_it = reg.finite_vertices_begin();
         dv_it != reg.finite_vertices_end();
         ++dv_it)
    {
      Reg_Vertex_handle dv = dv_it;
      os << dv->point().point()
         << " " << dv->point().weight()
         << " " << &*dv
         << std::endl;
    }
    os << std::endl;

    // write faces
    os << reg.number_of_faces() << std::endl;
    for (Reg_Finite_faces_iterator df_it = reg.finite_faces_begin();
         df_it != reg.finite_faces_end();
         ++df_it)
    {
      Reg_Face_handle df = df_it;
      for (int i=0; i<3; i++) {
        os << df->vertex(i)->point().point() << " ";
      }
      os << &*df
         << " " << &*df->neighbor(0) << "_" << &*df->neighbor(1) << "_" << &*df->neighbor(2)
         << " " << &*df->vertex(0) << "_" << &*df->vertex(1) << "_" << &*df->vertex(2)
         << std::endl;
    }
    os << std::endl;

    // we do not write neighbors

    os.close();
    if (verbose) {
      std::cerr << "INFO: Regular triangulation written to file \"" << filename << "\"" << std::endl;
    }
  } // >>>

// unit tests:
# define XSTR(s) STR(s)
# define STR(s) #s
# define TEST(x) if (not (x)) {fprintf(stderr, "assertion \"%s\" failed: file \"%s\", line %d\n", #x, __FILE__, __LINE__); abort();}
  void unit_tests() // <<<
  {
    unit_test_vor_faces_contain_duals();
    unit_test_vor_face_contains_point();
    unit_test_vor_locate();
    unit_test_vor_find_all_clusters();
    unit_test_vor_locate_in_clusters();
    unit_test_vor_cavity_measures();
    unit_test_treat_available_volumes_after_takeout();
    unit_test_treat_available_volumes();
  } // >>>
  void unit_test_vor_faces_contain_duals() // <<<
  // according to notes p.3276, the limit here is extra_radius <= 1.5, which is confirmed (1.49 OK, 1.51 NOT OK)
  {
    read_points("data/dump.gz");
    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 1.49);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 1.51);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(not vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
    { int n = 20;
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 1.49);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 1.51);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(not vor_faces_contain_duals(tri));
      TEST(vor_faces_are_convex(tri));
    }
  } // >>>
  void unit_test_vor_face_contains_point() // <<<
  // see notes p.3303 for a picture
  {
    read_points("data/dump.gz");
    Del_Vertex_handle dv;
    FT w = _domain.xmax() - _domain.xmin();
    FT h = _domain.ymax() - _domain.ymin();
    FT px(-3.4), py(-2.9); // point to find
    FT dvx(-3.02551), dvy(2.25544); // corresponding dv
    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      dv = del_vertex_close_to_point(Point(dvx, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(vor_face_contains_point(Point(px, py), dv, tri.get_offset(dv)+Shift(0,-1), tri));
      // just for fun: shift both point and dv
      TEST(vor_face_contains_point(Point(px+w, py+h), dv, tri.get_offset(dv)+Shift(1,0), tri));
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      dv = del_vertex_close_to_point(Point(dvx+w, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(vor_face_contains_point(Point(px+w, py+h), dv, tri.get_offset(dv)+Shift(0,0), tri));
      dv = del_vertex_close_to_point(Point(dvx, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(vor_face_contains_point(Point(px+3*w, py+3*h), dv, tri.get_offset(dv)+Shift(3,2), tri));
      dv = del_vertex_close_to_point(Point(dvx+2*w, dvy+2*h), tri); TEST(dv != Del_Vertex_handle());
      TEST(vor_face_contains_point(Point(px, py), dv, tri.get_offset(dv)+Shift(-2,-3), tri));
    }
  } // >>>
  void unit_test_vor_locate() // <<<
  // see notes p.3303 for a picture
  {
    std::pair<Del_Vertex_handle, Shift> res;
    read_points("data/dump.gz");
    Del_Vertex_handle dv;
    FT w = _domain.xmax() - _domain.xmin();
    FT h = _domain.ymax() - _domain.ymin();
    FT px(-3.4), py(-2.9); // point to find
    FT dvx(-3.02551), dvy(2.25544); // corresponding dv
    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      // find a point in the fundamental domain
      res = vor_locate(Point(px, py), tri);
      dv = del_vertex_close_to_point(Point(dvx, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(res.first == dv and res.second == Shift(0,-1));

      // find a point outside the fundamental domain
      res = vor_locate(Point(px+w, py+h), tri);
      dv = del_vertex_close_to_point(Point(dvx, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(res.first == dv and res.second == Shift(1,0));

      // find again a DelVertex
      TEST(vor_faces_contain_duals(tri));
      dv = del_vertex_close_to_point(Point(-2.67929,-2.19546), tri); TEST(dv != Del_Vertex_handle());
      res = vor_locate(dv->point(), tri);
      TEST(res.first == dv and res.second == Shift(0,0));
    }
    { int n = 20;
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      // make sure we are not kept in an infinite loop here...
      res = vor_locate(_orig_points[n], tri);
      TEST(squared_distance(_orig_points[n], unwrapped(res.first->point(), res.second)) > res.first->info().weight);
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      //del_save("delaunay.mytri", tri);
      //vor_save("voronoi.mytri", tri);

      // find a point in the fundamental domain
      res = vor_locate(Point(px, py), tri);
      dv = del_vertex_close_to_point(Point(dvx, dvy+2*h), tri); TEST(dv != Del_Vertex_handle());
      TEST(res.first == dv and res.second == tri.get_offset(dv)+Shift(0,-3));

      // find a point outside the fundamental domain
      res = vor_locate(Point(px+w, py+h), tri);
      dv = del_vertex_close_to_point(Point(dvx+w, dvy), tri); TEST(dv != Del_Vertex_handle());
      TEST(res.first == dv and res.second == tri.get_offset(dv)+Shift(0,0));

      // find again a DelVertex
      TEST(vor_faces_contain_duals(tri));
      dv = del_vertex_close_to_point(Point(-2.67929,-2.19546), tri); TEST(dv != Del_Vertex_handle());
      res = vor_locate(dv->point(), tri);
      TEST(res.first == dv and res.second == Shift(0,0));
    }
    { int n = 20;
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      // make sure we are not kept in an infinite loop here...
      res = vor_locate(_orig_points[n], tri);
      TEST(squared_distance(_orig_points[n], unwrapped(res.first->point(), res.second)) > res.first->info().weight);
    }
  } // >>>
  void unit_test_vor_find_all_clusters() // <<<
  {
    read_points("data/dump.gz");
    Cluster cl;
    Clusters clusters;
    DF_TO_CLS df_to_cls;
    //FT missed;
    Point point;
    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == tri.number_of_faces());
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 0.06);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 35 and clusters.size() == tri.number_of_faces());
    }
    { int n = 4; // 4 or 0
      // see what happens for a radius extension?
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      vor_find_all_clusters(clusters, df_to_cls, tri);
      //del_save("delaunay.mytri", tri);
      //vor_save("voronoi.mytri", tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == (n==4?3:1));
    }
    { // see what happens for a really large extension
      Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 5.0);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 0 and clusters.size() == 0);
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == tri.number_of_faces());
    }
    { Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 0.06);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 35 and clusters.size() == tri.number_of_faces());
    }
    { int n = 4; // 4 or 0
      // see what happens for a radius extension?
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == (n==4?3:1));
    }
    { // see what happens for a really large extension
      Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 5.0);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));

      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 0 and clusters.size() == 0);
    }
  } // >>>
  void unit_test_vor_locate_in_clusters() // <<<
  // For a picture, see notes p.3306
  {
    read_points("data/dump.gz");
    FT w = _domain.xmax() - _domain.xmin();
    FT h = _domain.ymax() - _domain.ymin();
    Clusters clusters;
    DF_TO_CLS df_to_cls;
    Del_Vertex_handle dv;
    Del_Face_handle df;
    Cluster cl;
    Point p;
    //FT missed;
    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 0.06);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 35 and clusters.size() == tri.number_of_faces());
      //del_save("delaunay.mytri", tri);
      //vor_save("voronoi.mytri", tri);

      // Delaunay points are not in any cluster
      dv = del_vertex_close_to_point(Point(-2.67929,-2.19546), tri); TEST(dv != Del_Vertex_handle());
      cl = vor_locate_in_clusters(dv->point(), clusters, df_to_cls, tri);
      TEST(cl == Cluster(NULL));

      // Find again a Voronoi vertex itself
      p = Point(-3.21894, -2.6065); // one of the Voronoi vertices
      df = vor_vertex_close_to_point(Point(p.x()+w, p.y()+h), tri); TEST(df != Del_Face_handle());
      //cl = vor_locate_in_clusters(df->info().point, clusters, df_to_cls, tri);
      cl = vor_locate_in_clusters(Point(p.x()-0.00, p.y()-0.01), clusters, df_to_cls, tri);
      TEST(cl == clusters.find(df_to_cls[df]));
    }
    { int n = 1;
      // see what happens for a radius extension:
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() >= 1);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl != Cluster(NULL));
    }
    { int n = 2;
      // see what happens for a large extension:
      Tri tri = del_from_nonperiodic_regular(n, 15.0);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_not_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 0 and clusters.size() == 0);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl == Cluster(NULL));
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular(SIZE_T_NONE, 0.06);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 35 and clusters.size() == tri.number_of_faces());
      //del_save("delaunay.mytri", tri);
      //vor_save("voronoi.mytri", tri);

      // Delaunay points are not in any cluster
      dv = del_vertex_close_to_point(Point(-2.67929,-2.19546), tri); TEST(dv != Del_Vertex_handle());
      cl = vor_locate_in_clusters(dv->point(), clusters, df_to_cls, tri);
      TEST(cl == Cluster(NULL));

      // Find again a Voronoi vertex itself
      p = Point(-3.21894, -2.6065); // one of the Voronoi vertices
      df = vor_vertex_close_to_point(Point(p.x()+3*w, p.y()+3*h), tri); TEST(df != Del_Face_handle());
      //cl = vor_locate_in_clusters(df->info().point, clusters, df_to_cls, tri);
      cl = vor_locate_in_clusters(Point(p.x()-0.00, p.y()-0.01), clusters, df_to_cls, tri);
      TEST(cl == clusters.find(df_to_cls[df]));

      // not allowed to search for points outside the fundamental domain
      cl = vor_locate_in_clusters(Point(p.x()+w, p.y()+h), clusters, df_to_cls, tri);
      TEST(cl == clusters.find(df_to_cls[df]));
    }
    { int n = 1;
      // see what happens for a radius extension:
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() >= 1);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl != Cluster(NULL));
    }
    { int n = 2;
      // see what happens for a large extension:
      Tri tri = del_from_nonperiodic_regular(n, 15.0);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_not_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 0 and clusters.size() == 0);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl == Cluster(NULL));
    }
  } // >>>
  void unit_test_vor_cavity_measures() // <<<
  {
    read_points("data/dump.gz");
    FT w = _domain.xmax() - _domain.xmin();
    FT h = _domain.ymax() - _domain.ymin();
    Clusters clusters;
    Cluster cl;
    DF_TO_CLS df_to_cls;
    Cavity_measures cms;
    Cavity_measure cm;
    FT eps(1.0e-5);

    { // 1-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == tri.number_of_faces());

      cms = vor_cavity_measures(tri, clusters);
      TEST(cms.size() == 1);
      cm = cms.begin()->second;
      TEST(CGAL::abs(cm.vol_triangles - w*h) < eps);
      TEST(CGAL::abs(cm.vol_cavity - (w*h - del_vertices_spheresvolume(tri))) < eps);
      TEST(CGAL::abs(cm.surf - del_vertices_spheressurface(tri)) < eps);
    }
    { int n = 0;
      // see what happens for a radius extension:
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_1_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == 1);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl != Cluster(NULL));

      cms = vor_cavity_measures(tri, clusters, cl);
      TEST(cms.size() == 1);
      cm = cms.begin()->second;
      TEST(cm.vol_triangles > cm.vol_cavity);
      TEST(cm.surf > FT(0));
    }

    { // 9-sheeted
      Tri tri = del_from_nonperiodic_regular();
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      TEST(vor_vertices_number_of_covered(tri) == 0);
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() == tri.number_of_faces());
      //del_save("delaunay.mytri", tri);
      //vor_save("voronoi.mytri", tri);

      cms = vor_cavity_measures(tri, clusters);
      TEST(cms.size() == 1);
      cm = cms.begin()->second;
      TEST(CGAL::abs(cm.vol_triangles - w*h) < eps);
      TEST(CGAL::abs(cm.vol_cavity - (w*h - del_vertices_spheresvolume(tri))) < eps);
      TEST(CGAL::abs(cm.surf - del_vertices_spheressurface(tri)) < eps);
    }
    { int n = 1;
      // see what happens for a radius extension:
      Tri tri = del_from_nonperiodic_regular(n, _orig_radii[n]);
      tri.convert_to_9_sheeted_covering();
      vor_vertices_create(tri);
      TEST(vor_faces_are_convex(tri));
      vor_find_all_clusters(clusters, df_to_cls, tri);
      TEST(clusters.number_of_sets() == 1 and clusters.size() >= 1);

      cl = vor_locate_in_clusters(_orig_points[n], clusters, df_to_cls, tri);
      TEST(cl != Cluster(NULL));

      cms = vor_cavity_measures(tri, clusters, cl);
      TEST(cms.size() == 1);
      cm = cms.begin()->second;
      TEST(cm.vol_triangles > cm.vol_cavity);
      TEST(cm.surf > FT(0));
    }
  } // >>>
  void unit_test_treat_available_volumes_after_takeout() // <<<
  {
    read_points("data/dump.gz");
    std::cout << "Free volumes in data/dump.gz" << std::endl;
    // We use only one cavity because it tests more functions
    treat_available_volumes_after_takeout(true);
  } // >>>
  void unit_test_treat_available_volumes() // <<<
  {
    read_points("data/dump.gz");
    std::cout << "Available volumes in data/dump.gz" << std::endl;
    // We use only one cavity because it tests more functions
    treat_available_volumes();
  } // >>>

// convenience functions for unit tests:
  Del_Vertex_handle del_vertex_close_to_point(const Point& p, const Tri& tri, const FT& dist=1.0e-2) const // <<<
  {
    for (Del_Vertex_iterator dv_it = tri.vertices_begin();
         dv_it != tri.vertices_end();
         ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      if (squared_distance(unwrapped(dv->point(), tri.get_offset(dv)), p) < dist*dist) {
        return dv;
      }
    }
    return Del_Vertex_handle();
  } // >>>
  Del_Face_handle vor_vertex_close_to_point(const Point& p, const Tri& tri, const FT& dist=1.0e-2) const // <<<
  {
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it;
      if (squared_distance(unwrapped(df->info().point, df->info().offset), p) < dist*dist) {
        return df;
      }
    }
    return Del_Face_handle();
  } // >>>
  FT del_vertices_spheresvolume(const Tri& tri) const // <<<
  {
    FT result(0);
    for (Del_Unique_vertex_iterator dv_it = tri.unique_vertices_begin(),
                                    dv_end = tri.unique_vertices_end();
         dv_it != dv_end; ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      result += dv->info().weight;
    }
    return result*FT(M_PI);
  } // >>>
  FT del_vertices_spheressurface(const Tri& tri) const // <<<
  {
    FT result(0);
    for (Del_Unique_vertex_iterator dv_it = tri.unique_vertices_begin(),
                                    dv_end = tri.unique_vertices_end();
         dv_it != dv_end; ++dv_it)
    {
      Del_Vertex_handle dv = dv_it;
      result += CGAL::sqrt(dv->info().weight);
    }
    return FT(2*M_PI)*result;
  } // >>>
  size_t vor_vertices_number_of_covered(const Tri& tri) const // <<<
  {
    size_t result = 0;
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it; // cast operator
      const Vor_vertex& vv = df->info();
      if (vv.dual_is_canonical and vv.is_covered)
        ++result;
    }
    return result;
  } // >>>
  size_t vor_vertices_number_of_not_covered(const Tri& tri) const // <<<
  {
    size_t result = 0;
    for (Del_Face_iterator df_it = tri.faces_begin();
         df_it != tri.faces_end();
         ++df_it)
    {
      Del_Face_handle df = df_it; // cast operator
      const Vor_vertex& vv = df->info();
      if (vv.dual_is_canonical and not vv.is_covered)
        ++result;
    }
    return result;
  } // >>>

};

# endif // PERIODIC_2_DUAL_TRIANGULATIONS_2_H
// vim:foldmethod=marker:foldmarker=<<<,>>>
