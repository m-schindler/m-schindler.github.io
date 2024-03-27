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
# ifndef TOOLS_H
# define TOOLS_H

template<class T>
T signum(const T t) // <<<
{
  if (t == 0) return T(0);
  else return ((t < 0) ? T(-1) : T(1));
} // >>>

struct Ptr_hash_function { // <<<
  // uses t.ptr()
  typedef std::size_t result_type;
  template<class T>
  std::size_t operator() (const T& t) const {
    return std::size_t(t.ptr());
  }
}; // >>>

struct Handle_hash_function { // <<<
  // drop-in replacement for CGAL::Handle_hash_function
  // which did not work reliably in del_unique_face_hash
  typedef std::size_t result_type;
  template<class T>
  std::size_t operator() (const T& t) const {
    return std::size_t(&*t);
  }
}; // >>>

template<class Hhash>
struct Handle_Offset_2_hash_function { // <<<
  // hash function for the pair<Handle, Offset>
  typedef size_t result_type;
  Hhash handle_hash;
  std::hash<int> int_hash;

  template<class Handle, class Offset>
  size_t operator() (const std::pair<Handle, Offset>& p) const
  {
    // combine the hash values:
    // this is the implementation of boost::combine_hash
    // CGAL::Handle_hash_function divides the memory address by the sizeof
    // XXX is it consistent to combine CGAL::Handle_hash_function with std::hash<int>?
    //     At least, I had problems with this combination
    size_t
    result  = handle_hash(p.first)   + 0x9e3779b9;
    result ^= int_hash(p.second.x()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(p.second.y()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    return result;
  }
}; // >>>

template<class Hhash>
struct Handle_Offset_3_hash_function { // <<<
  // hash function for the pair<Handle, Offset>
  typedef size_t result_type;
  Hhash handle_hash;
  std::hash<int> int_hash;

  template<class Handle, class Offset>
  size_t operator() (const std::pair<Handle, Offset>& p) const
  {
    // combine the hash values:
    // this is the implementation of boost::combine_hash
    // CGAL::Handle_hash_function divides the memory address by the sizeof
    // XXX is it consistent to combine CGAL::Handle_hash_function with std::hash<int>?
    //     At least, I had problems with this combination
    size_t
    result  = handle_hash(p.first)   + 0x9e3779b9;
    result ^= int_hash(p.second.x()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(p.second.y()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    result ^= int_hash(p.second.z()) + 0x9e3779b9 + (result << 6) + (result >> 2);
    return result;
  }
}; // >>>

template<class FT>
class Cavity_measure_base // <<<
{
  public:
  FT vol_triangles, vol_cavity, surf, surf_radii;
  Cavity_measure_base() : vol_triangles(0), vol_cavity(0), surf(0), surf_radii(0) {}
  void clear()
  {
    vol_triangles = 0;
    vol_cavity = 0;
    surf = 0;
    surf_radii = 0;
  }
  void add(const FT vol_total, const FT vol_cav, const FT srf, const FT radii)
  {
    vol_triangles += vol_total;
    vol_cavity += vol_cav;
    surf += srf;
    surf_radii += srf*radii;
  }
}; // >>>

template<class vertid_t, class shift_t>
inline bool tuple_less(const vertid_t this_v, const shift_t this_o[NDIM], const vertid_t other_v, const shift_t other_o[NDIM]) // <<<
// general tuple comparison on (vid, o0, o1, o2)
//   first entry dominates, then second, ...
//   as if the first were multiplied with a large number, the second with a smaller, the third with an even smaller, ...
//   such that the range of each possible value in each entry is outrun by all the factors.
{
  if (this_v < other_v)
    return true;
  else if (this_v > other_v)
    return false;

  for (unsigned int i=0; i<NDIM; i++)
    if (this_o[i] < other_o[i])
      return true;
    else if (this_o[i] > other_o[i])
      return false;

  return false; // they are equal
} // >>>
template<class VERTID_t, class SHIFT_t>
class Periodic_simplex_fingerprint // <<<
// a hash-like object which allows to uniquely oder periodic simplices
// (used in std::map)
{
private:
  VERTID_t verts[SN];
  SHIFT_t shs[NDIM][NDIM];
public:
  typedef VERTID_t vertid_t;
  typedef SHIFT_t shift_t;
  // constructor:
  Periodic_simplex_fingerprint(const vertid_t v[SN], const shift_t o[SN][NDIM])
  {
    # if SN != NDIM + 1
    #   error "Incompatible SN and NDIM"
    # endif

    // sort vertices and shifts:
    # if NDIM == 3
    unsigned int i=0, j=1, k=2, l=3;
    if (tuple_less(v[l],o[l], v[k],o[k])) std::swap(l,k);
    if (tuple_less(v[k],o[k], v[j],o[j])) std::swap(k,j);
    if (tuple_less(v[l],o[l], v[k],o[k])) std::swap(l,k);
    if (tuple_less(v[j],o[j], v[i],o[i])) std::swap(j,i);
    if (tuple_less(v[k],o[k], v[j],o[j])) std::swap(k,j);
    if (tuple_less(v[l],o[l], v[k],o[k])) std::swap(l,k);
    // must be ordered now (and they may not be equal)
    //assert(not tuple_less(v[j],o[j], v[i],o[i]));
    //assert(not tuple_less(v[k],o[k], v[j],o[j]));
    //assert(not tuple_less(v[l],o[l], v[k],o[k]));
    //// they may not be equal:
    //if (not (tuple_less(v[i],o[i], v[j],o[j]) and tuple_less(v[j],o[j], v[k],o[k]) and tuple_less(v[k],o[k], v[l],o[l]))) {
    //  std::cerr << "i=" << i << " v[i]=" << v[i] << " o[i]=(" << o[i][0] << " " << o[i][1] << " " << o[i][2] << ")" << std::endl;
    //  std::cerr << "j=" << j << " v[j]=" << v[j] << " o[j]=(" << o[j][0] << " " << o[j][1] << " " << o[j][2] << ")" << std::endl;
    //  std::cerr << "k=" << k << " v[k]=" << v[k] << " o[k]=(" << o[k][0] << " " << o[k][1] << " " << o[k][2] << ")" << std::endl;
    //  std::cerr << "l=" << l << " v[l]=" << v[l] << " o[l]=(" << o[l][0] << " " << o[l][1] << " " << o[l][2] << ")" << std::endl;
    //}
    assert(tuple_less(v[i],o[i], v[j],o[j]));
    assert(tuple_less(v[j],o[j], v[k],o[k]));
    assert(tuple_less(v[k],o[k], v[l],o[l]));
    # elif NDIM == 2
    unsigned int i=0, j=1, k=2;
    if (tuple_less(v[k],o[k], v[j],o[j])) std::swap(k,j);
    if (tuple_less(v[j],o[j], v[i],o[i])) std::swap(j,i);
    if (tuple_less(v[k],o[k], v[j],o[j])) std::swap(k,j);
    assert(tuple_less(v[i],o[i], v[j],o[j]));
    assert(tuple_less(v[j],o[j], v[k],o[k]));
    # endif

    // vertices
    verts[0] = v[i];
    verts[1] = v[j];
    verts[2] = v[k];
    # if NDIM >= 3
    verts[3] = v[l];
    # endif

    // relative offsets (shifts)
    for (unsigned int dim=0; dim<NDIM; dim++) {
      shs[0][dim] = o[j][dim] - o[i][dim];
      shs[1][dim] = o[k][dim] - o[i][dim];
      # if NDIM >= 3
      shs[2][dim] = o[l][dim] - o[i][dim];
      # endif
    }
  }

  bool operator<(const Periodic_simplex_fingerprint& other) const
  // general tuple comparison on the tuple (v[0],v[1],v[2],v[3], o[0][0],o[0][1],...,o[2][1],o[2][2])
  //   first entry dominates, then second, ...
  //   as if the first were multiplied with a large number, the second with a smaller, the third with an even smaller, ...
  //   such that the range of each possible value in each entry is outrun by the factors.
  {
    for (unsigned int n=0; n<SN; n++) {
      if (this->verts[n] < other.verts[n]) {
        return true;
      } else if (this->verts[n] > other.verts[n]) {
        return false;
      }
    }

    for (unsigned int n=0; n<NDIM; n++) {
      for (unsigned int dim=0; dim<NDIM; dim++) {
        if (this->shs[n][dim] < other.shs[n][dim]) {
          return true;
        } else if (this->shs[n][dim] > other.shs[n][dim]) {
          return false;
        }
      }
    }

    // they are equal
    return false;
  }

  // TODO: include stream correctly
  //std::ostream& operator<<(std::ostream& stream) const
  void print(std::ostream& stream) const
  {
    stream << "fingerprint(";
    for (unsigned int n=0; n<SN; n++)
      stream << verts[n] << ",";
    stream << " ";
    for (unsigned int n=0; n<NDIM; n++) {
      for (unsigned int dim=0; dim<NDIM; dim++) {
        stream << shs[n][dim] << ",";
      }
      stream << " ";
    }
    stream << ")";
    stream << std::endl;
    //return stream;
  }

}; // >>>

template<class Vertex_handle, class Offset, class Point>
class Periodic_vertex // <<<
// pair of a fundamental vertex and an offset
{
public:
  Periodic_vertex(const Vertex_handle& fv, const Offset& o)
    : fund_vh(fv), off(o)
    {}

  const Vertex_handle& vertex_handle() const {return fund_vh;}
  const Offset& offset() const {return off;}

  bool operator<(const Periodic_vertex& other) const
  // general tuple comparison in the order x, y, z:
  // for each coordinate, compare the tuple (offset, point) because the points are within the original domain
  // TODO: for some usage, it would be enough to compare memory addresses of vertices
  //       is this identical in all used cases?
  {
    # if 0
    Point thispt = this->fund_vh->point(),
          otherpt = other.fund_vh->point();

    for (unsigned int dim=0; dim<NDIM; dim++) {
      // test offsets:
      if (this->off[dim] < other.off[dim]) return true;
      else if (this->off[dim] > other.off[dim]) return false;
      // the offsets are equal in this dimension

      // test vertices:
      if (this->fund_vh != other.fund_vh) {
        if (thispt[dim] < otherpt[dim]) return true;
        else if (thispt[dim] > otherpt[dim]) return false;
      }
      // the vertices are equal
    }
    # else
    for (unsigned int dim=0; dim<NDIM; dim++) {
      // test offsets:
      if (this->off[dim] < other.off[dim]) return true;
      else if (this->off[dim] > other.off[dim]) return false;
    }
    // if we get here, the offsets are equal

    // test memory address (hash) of the vertices themselves:
    // TODO use hashes provided by STL_extension (how to get that from git?)
    const size_t this_hash = size_t(&*(this->fund_vh)),
                other_hash = size_t(&*(other.fund_vh));
    if (this_hash < other_hash) return true;
    else if (this_hash > other_hash) return false;
    # endif

    // they are equal
    return false;
  }

  bool operator==(const Periodic_vertex& other) const
  // general tuple comparison (off_x, off_y, off_z, fund_vh)
  {
    for (unsigned int dim=0; dim<NDIM; dim++)
      if (this->off[dim] != other.off[dim]) return false;

    // test memory address (hash) of the vertices themselves:
    // TODO use hashes provided by STL_extension (how to get that from git?)
    const size_t this_hash = size_t(&*(this->fund_vh)),
                other_hash = size_t(&*(other.fund_vh));
    return (this_hash == other_hash);
  }

private:
  Vertex_handle fund_vh;
  Offset off;
}; // >>>

template<class Cell_handle, class Shift>
class Periodic_cell // <<<
// pair of a canonical cell and an offset
{
public:
  // default constructor to allow empty cell_handle
  Periodic_cell()
    : canon_ch(Cell_handle()), sh(Shift()) {}

  Periodic_cell(const Cell_handle& cch, const Shift& shift)
    : canon_ch(cch), sh(shift)
    {
      assert(canon_ch->info().dual_is_canonical);
    }

  const Cell_handle& cell_handle() const {return canon_ch;}
  const Shift& shift() const {return sh;}

  bool operator<(const Periodic_cell& other) const
  // general tuple comparison (sh_x, sh_y, sh_z, canon_ch)
  {
    for (unsigned int dim=0; dim<NDIM; dim++) {
      if (this->sh[dim] < other.sh[dim]) return true;
      else if (this->sh[dim] > other.sh[dim]) return false;
    }
    // if we get here, the shifts are equal

    // test memory address (hash) of the cells themselves:
    // TODO use hashes provided by STL_extension (how to get that from git?)
    const size_t this_hash = size_t(&*(this->canon_ch)),
                other_hash = size_t(&*(other.canon_ch));
    if (this_hash < other_hash) return true;
    else if (this_hash > other_hash) return false;

    // if we get here, also the cells are equal
    return false;
  }

  bool operator==(const Periodic_cell& other) const
  // general tuple comparison (sh_x, sh_y, sh_z, canon_ch)
  {
    for (unsigned int dim=0; dim<NDIM; dim++)
      if (this->sh[dim] != other.sh[dim]) return false;

    // test memory address (hash) of the cells themselves:
    // TODO use hashes provided by STL_extension (how to get that from git?)
    const size_t this_hash = size_t(&*(this->canon_ch)),
                other_hash = size_t(&*(other.canon_ch));
    return (this_hash == other_hash);
  }

private:
  Cell_handle canon_ch;
  Shift sh;
}; // >>>

# endif // TOOLS_H
// vim:foldmethod=marker:foldmarker=<<<,>>>
