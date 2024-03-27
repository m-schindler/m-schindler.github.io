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
# ifndef PERIODIC_TRIANGULATION_COVER_ITERATORS_2_H
# define PERIODIC_TRIANGULATION_COVER_ITERATORS_2_H
# include <list>

template<class PT>
class Periodic_face_cover_iterator_2
// similar to Triangle_iterator(COVER)
// This iterator runs over all faces that cover a part of the
// domain (1-sheeted or 9-sheeted)
{
  typedef typename PT::Face_handle       Face_handle;
  typedef typename PT::Face_iterator     Face_iterator;
  typedef typename PT::Offset            Offset;
  typedef typename PT::Offset            Shift;

  const PT&                     pt;
  Shift sh0, shx, shy, shxy;
  std::list<Shift*>              shifts;
  typename std::list<Shift*>::iterator    s_it;
  Face_iterator                 f_it;

  void reset_shifts(Face_handle f)
  {
    shifts.clear();

    // add new shifts:
    Offset off0 = pt.int_to_off(f->offset(0)),
           off1 = pt.int_to_off(f->offset(1)),
           off2 = pt.int_to_off(f->offset(2));
    assert(off0.x()==0 or off1.x()==0 or off2.x()==0);
    assert(off0.y()==0 or off1.y()==0 or off2.y()==0);
    bool cross_x = (off0.x() != off1.x() or off0.x() != off2.x() or off1.x() != off2.x());
    bool cross_y = (off0.y() != off1.y() or off0.y() != off2.y() or off1.y() != off2.y());
    if (cross_x) {
      shifts.insert(shifts.begin(), &shx);
      if (cross_y) shifts.insert(shifts.begin(), &shxy);
    }
    if (cross_y) shifts.insert(shifts.begin(), &shy);
    shifts.insert(shifts.begin(), &sh0);
  }

  public:
  Periodic_face_cover_iterator_2(const PT& pt_)
  : pt(pt_),
    sh0(Shift(0, 0)),
    shx(Shift(-pt_.number_of_sheets()[0], 0)),
    shy(Shift(0, -pt_.number_of_sheets()[1])),
    shxy(Shift(-pt_.number_of_sheets()[0], -pt_.number_of_sheets()[1]))
  {
    f_it = pt.faces_begin();
    Face_handle f = f_it; //cast
    reset_shifts(f);
    s_it = shifts.begin();
  }

  Periodic_face_cover_iterator_2& operator++()
  {
    ++s_it;
    if (s_it == shifts.end()) {
      ++f_it;
      if (f_it != pt.faces_end()) {
        Face_handle f = f_it; //cast
        reset_shifts(f);
        s_it = shifts.begin();
      }
    }
    return *this;
  }

  Face_handle face_handle() const
  {
    Face_handle f = f_it;
    return f;
  }
  Shift* face_shift() const {return *s_it;}

  // cast operators
  //operator const Face_handle () const {return *f_it;}
  //operator const Shift () const {return *s_it;}

  bool atend() const
  {
    return f_it == pt.faces_end() and s_it == shifts.end();
  }
};


# endif // PERIODIC_TRIANGULATION_COVER_ITERATORS_2_H
