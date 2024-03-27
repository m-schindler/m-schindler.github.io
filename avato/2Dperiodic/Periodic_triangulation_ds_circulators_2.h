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
// This file is takes up and modifies some code from CGAL/Triangulation_data_structure_2.h.
#ifndef PERIODIC_TRIANGULATION_DS_CIRCULATORS_2_H
#define PERIODIC_TRIANGULATION_DS_CIRCULATORS_2_H

#include <utility>
#include <iterator>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

template < class Tds>
class Periodic_triangulation_ds_face_circulator_2
  : public CGAL::Bidirectional_circulator_base< typename Tds::Face,
                                 std::ptrdiff_t,
                                 std::size_t>,
    public CGAL::Triangulation_cw_ccw_2

{
// Modifies the face_circulator such that it works even in
// periodic triangulations which are 1-sheeted but not "valid ones":
// A vertex may occur several times in a face, with different offsets
// -- such a circulator is required when constructing a periodic voronoi
// "triangulation" from another (Delaunay/Regular) one
private:
  typedef
  CGAL::Bidirectional_circulator_base< typename Tds::Face,
                                 std::ptrdiff_t,
                                 std::size_t>  Base_circulator;

public:
  typedef Periodic_triangulation_ds_face_circulator_2<Tds> Face_circulator;
  typedef typename Tds::Face                      Face;
  typedef typename Tds::Vertex                    Vertex;
  typedef typename Tds::Face_handle               Face_handle;
  typedef typename Tds::Vertex_handle             Vertex_handle;


private:
  Vertex_handle _v;
  Face_handle    pos;
  int ind; // index of center vertex in the face handle

public:
  Periodic_triangulation_ds_face_circulator_2()
    : _v(), pos(), ind(-1)
  {}

  Periodic_triangulation_ds_face_circulator_2(Vertex_handle v,
                                     Face_handle f = Face_handle(),
                                     int i = -1);

  Face_circulator& operator=(const Face_circulator& other);
  Face_circulator& operator++();
  Face_circulator operator++(int);
  Face_circulator& operator--();
  Face_circulator operator--(int);

  bool operator==(const Face_circulator &fc) const;
  bool operator!=(const Face_circulator &fc) const;

  bool is_empty() const;
  // what are these two for?
  //bool operator==(Nullptr_t CGAL_triangulation_assertion_code(n)) const;
  //bool operator!=(Nullptr_t CGAL_triangulation_assertion_code(n)) const;

  Face_handle base()  const {return pos;}
  operator Face_handle()  const {return pos;} // cast operator
  int index() const {return ind;}
  operator int()  const {return ind;} // cast operator
};

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds> ::
Periodic_triangulation_ds_face_circulator_2(Vertex_handle v, Face_handle f, int i)
  : _v(v), pos(f), ind(i)
{
  if (_v == Vertex_handle()) { pos = Face_handle(); ind = -1; }
  else if ( pos == Face_handle()) {
    pos = v->face();
    if (ind == -1) ind = pos->index(v); // choose *one* occurrence of v in pos
    else if (pos->vertex(ind) != _v) ind = pos->index(v); // choose *one* occurrence of v in pos
  }

  if (pos == Face_handle() || pos->dimension() < 2) {
    _v = Vertex_handle() ; pos = Face_handle(); ind = -1; return;}
  else CGAL_triangulation_precondition(pos->has_vertex(v) and pos->vertex(ind) == v);
}

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds>&
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator=(const Periodic_triangulation_ds_face_circulator_2<Tds>& other)
{
   _v = other._v;
  pos = other.pos;
  ind = other.ind;
  return *this;
}

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds>&
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator++()
{
  CGAL_triangulation_precondition(pos != Face_handle() &&
                                 _v != Vertex_handle() &&
                                 ind != -1 && pos->vertex(ind) == _v);
  // get the counter-clockwise nbr of pos around the _v
  // turn ccw on the vertices of f and then determine nbr
  Face_handle nbr = pos->neighbor(ccw(ind));
  // identify the index of the very vertex, starting with the nbr
  int nbrind = ccw(nbr->index(pos));
  CGAL_expensive_assertion(nbr->neighbor(nbrind) != pos and nbr->neighbor(ccw(nbrind)) != pos);
  CGAL_expensive_assertion(nbr->vertex(nbrind) == _v);
  ind = nbrind;
  pos = nbr;

  return *this;
}

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds>
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator++(int)
{
  Face_circulator tmp(*this);
  ++(*this);
  return tmp;
}

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds>&
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator--()
{
  CGAL_triangulation_precondition(pos != Face_handle() &&
                                 _v != Vertex_handle() &&
                                 ind != -1 && pos->vertex(ind) == _v);
  Face_handle nbr = pos->neighbor(cw(ind));
  int nbrind = cw(nbr->index(pos));
  CGAL_expensive_assertion(nbr->neighbor(nbrind) != pos and nbr->neighbor(cw(nbrind)) != pos);
  CGAL_expensive_assertion(nbr->vertex(nbrind) == _v);
  ind = nbrind;
  pos = nbr;

  return *this;
}

template < class Tds >
Periodic_triangulation_ds_face_circulator_2<Tds>
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator--(int)
{
  Face_circulator tmp(*this);
  --(*this);
  return tmp;
}

template < class Tds >
inline bool
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator==(const Periodic_triangulation_ds_face_circulator_2<Tds>&fc) const
{
  return (_v == fc._v) &&  (pos == fc.pos) && (ind == fc.ind);
}

template < class Tds >
inline bool
Periodic_triangulation_ds_face_circulator_2<Tds> ::
operator!=(const Periodic_triangulation_ds_face_circulator_2<Tds>&fc) const
{
return ! (*this == fc);
}

template < class Tds >
inline bool
Periodic_triangulation_ds_face_circulator_2<Tds> ::
is_empty() const
{
return (_v == Vertex_handle() ||  pos == Face_handle() || ind == -1);
}

# endif // PERIODIC_TRIANGULATION_DS_CIRCULATORS_2_H
