// Copyright (c) 2017  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Simon Giraudot

#ifndef CGAL_READ_LAS_POINTS_H
#define CGAL_READ_LAS_POINTS_H

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <lasreader_las.hpp>
#pragma GCC diagnostic pop
#endif

namespace CGAL {


namespace LAS
{

  /// \cond SKIP_IN_MANUAL
  namespace Property
  {
    
    struct X                   { typedef double type; };
    struct Y                   { typedef double type; };
    struct Z                   { typedef double type; };
    struct intensity           { typedef unsigned short type; };
    struct return_number       { typedef unsigned char type; };
    struct number_of_returns   { typedef unsigned char type; };
    struct scan_direction_flag { typedef unsigned char type; };
    struct edge_of_flight_line { typedef unsigned char type; };
    struct classification      { typedef unsigned char type; };
    struct synthetic_flag      { typedef unsigned char type; };
    struct keypoint_flag       { typedef unsigned char type; };
    struct withheld_flag       { typedef unsigned char type; };
    struct scan_angle          { typedef float type; };
    struct user_data           { typedef unsigned char type; };
    struct point_source_ID     { typedef unsigned short type; };
    struct deleted_flag        { typedef unsigned int type; };
    struct gps_time            { typedef double type; };
    struct R                   { typedef unsigned short type; };
    struct G                   { typedef unsigned short type; };
    struct B                   { typedef unsigned short type; };
    struct I                   { typedef unsigned short type; };
  }
  /// \endcond


  /**
     \ingroup PkgPointSetProcessing
     
     Generates a LAS property handler to read 3D points. Points are
     constructed from the input the using 3 LAS properties
     `LAS::Property::X`, `LAS::Property::Y` and `LAS::Property::Z`.

     \sa `read_las_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  cpp11::tuple<PointMap,
               typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
               Property::X, Property::Y, Property::Z >
  point_reader(PointMap point_map)
  {
    return cpp11::make_tuple (point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
                              Property::X(), Property::Y(), Property::Z());
  }
  
/// \cond SKIP_IN_MANUAL
  
namespace internal {

  void get_value(const LASpoint& r, double& v, LAS::Property::X&)
  { v = r.get_x(); }
  void get_value(const LASpoint& r, double& v, LAS::Property::Y&)
  { v = r.get_y(); }
  void get_value(const LASpoint& r, double& v, LAS::Property::Z&)
  { v = r.get_z(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::intensity&)
  { v = r.get_intensity(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::return_number&)
  { v = r.get_return_number(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::number_of_returns&)
  { v = r.get_number_of_returns(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::scan_direction_flag&)
  { v = r.get_scan_direction_flag(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::edge_of_flight_line&)
  { v = r.get_edge_of_flight_line(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::classification&)
  { v = r.get_classification(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::synthetic_flag&)
  { v = r.get_synthetic_flag(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::keypoint_flag&)
  { v = r.get_keypoint_flag(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::withheld_flag&)
  { v = r.get_withheld_flag(); }
  void get_value(const LASpoint& r, float& v, LAS::Property::scan_angle&)
  { v = r.get_scan_angle(); }
  void get_value(const LASpoint& r, unsigned char& v, LAS::Property::user_data&)
  { v = r.get_user_data(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::point_source_ID&)
  { v = r.get_point_source_ID(); }
  void get_value(const LASpoint& r, unsigned int& v, LAS::Property::deleted_flag&)
  { v = r.get_deleted_flag(); }
  void get_value(const LASpoint& r, double& v, LAS::Property::gps_time&)
  { v = r.get_gps_time(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::R&)
  { v = r.get_R(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::G&)
  { v = r.get_G(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::B&)
  { v = r.get_B(); }
  void get_value(const LASpoint& r, unsigned short& v, LAS::Property::I&)
  { v = r.get_I(); }

  
  template <std::size_t N>
  struct Filler
  {
    template <class Value_tuple, class LAS_property_tuple>
    static void fill(const LASpoint& r, Value_tuple& values, LAS_property_tuple wrappers)
    {
      get_value(r, std::get<N>(values), std::get<N+2>(wrappers));
      Filler<N-1>::fill(r, values, wrappers);
    }
  };

  template<int ...>
  struct seq { };

  template<int N, int ...S>
  struct gens : gens<N-1, N-1, S...> { };

  template<int ...S>
  struct gens<0, S...> {
    typedef seq<S...> type;
  };

  template<class ValueType, class Functor, class Tuple, int ...S>
  ValueType call_functor(Functor f, Tuple t, seq<S...>) {
    return f(std::get<S>(t) ...);
  }

  template <class ValueType, class Functor, typename ... T>
  ValueType call_functor(Functor f, std::tuple<T...>& t)
  {
    return call_functor<ValueType>(f, t, typename gens<sizeof...(T)>::type());
  }

  template<>
  struct Filler<0>
  {
    template <class Value_tuple, class LAS_property_tuple>
    static void fill(const LASpoint& r, Value_tuple& values, LAS_property_tuple wrappers)
    {
      get_value(r, std::get<0>(values), std::get<2>(wrappers));
    }
  };
  
  template <typename OutputValueType, typename PropertyMap, typename T>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, T>& current);

  template <typename OutputValueType, typename PropertyMap, typename T,
            typename NextPropertyBinder, typename ... PropertyMapBinders>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, T>& current,
                           NextPropertyBinder& next,
                           PropertyMapBinders&& ... properties);
  
  template <typename OutputValueType,
            typename PropertyMap,
            typename Constructor,
            typename ... T>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::tuple<PropertyMap, Constructor, T...>& current)
  {
    typedef typename PropertyMap::value_type PmapValueType;
    std::tuple<typename T::type...> values;
    Filler<sizeof...(T)-1>::fill(reader, values, current);
    PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
    put (std::get<0>(current), new_element, new_value);
  }
  
  template <typename OutputValueType,
            typename PropertyMap,
            typename Constructor,
            typename ... T,
            typename NextPropertyBinder,
            typename ... PropertyMapBinders>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::tuple<PropertyMap, Constructor, T...>& current,
                           NextPropertyBinder& next,
                           PropertyMapBinders&& ... properties)
  {
    typedef typename PropertyMap::value_type PmapValueType;
    std::tuple<typename T::type...> values;
    Filler<sizeof...(T)-1>::fill(reader, values, current);
    PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
    put (std::get<0>(current), new_element, new_value);
  
    process_properties (reader, new_element, next, properties...);
  }


  template <typename OutputValueType, typename PropertyMap, typename T>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, T>& current)
  {
    typename T::type new_value = typename T::type();
    get_value (reader, new_value, current.second);
    put (current.first, new_element, new_value);
  }

  template <typename OutputValueType, typename PropertyMap, typename T,
            typename NextPropertyBinder, typename ... PropertyMapBinders>
  void process_properties (const LASpoint& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, T>& current,
                           NextPropertyBinder& next,
                           PropertyMapBinders&& ... properties)
  {
    typename T::type new_value = typename T::type();
    get_value (reader, new_value, current.second);
    put (current.first, new_element, new_value);
    process_properties (reader, new_element, next, properties...);
  }
  
} // namespace internal
  

/// \endcond

}

//===================================================================================
/// \ingroup PkgPointSetProcessing

/// Reads user-selected points properties from a .las or .laz stream.
/// Potential additional properties are ignored.
///
/// Properties are handled through a variadic list of property
/// handlers. A property handle can either be:
///
///  - A `std::pair<PropertyMap, LASProperty >` if the user wants to
///  read a LAS property as a scalar value `LASProperty::type` (for
///  example, storing an `int` LAS property into an `int` variable).
///
///  - A `CGAL::cpp11::tuple<PropertyMap, Constructor,
///  LASProperty...>` if the user wants to use one or several LAS
///  properties to construct a complex object (for example, storing 4
///  `unsigned short` LAS properties into a %Color object that can for
///  example be a `CGAL::cpp11::array<unsigned short, 4>`). In that
///  case, the second element of the tuple should be a functor that
///  constructs the value type of `PropertyMap` from N objects of
///  of type `LASProperty::type`.
///
/// The LAS standard defines a fixed set of properties accessible
/// through the following tag classes:
///
///  - `LAS::Property::X` with type `double`
///  - `LAS::Property::Y` with type `double`
///  - `LAS::Property::Z` with type `double`
///  - `LAS::Property::intensity` with type `unsigned short`
///  - `LAS::Property::return_number` with type `unsigned char`
///  - `LAS::Property::number_of_returns` with type `unsigned char`
///  - `LAS::Property::number_of_returns` with type `unsigned char`
///  - `LAS::Property::scan_direction_flag` with type `unsigned char`
///  - `LAS::Property::edge_of_flight_line` with type `unsigned char`
///  - `LAS::Property::classification` with type `unsigned char`
///  - `LAS::Property::synthetic_flag` with type `unsigned char`
///  - `LAS::Property::keypoint_flag` with type `unsigned char`
///  - `LAS::Property::withheld_flag` with type `unsigned char`
///  - `LAS::Property::scan_angle` with type `double`
///  - `LAS::Property::user_data` with type `unsigned char`
///  - `LAS::Property::point_source_ID` with type `unsigned short`
///  - `LAS::Property::deleted_flag` with type `unsigned int`
///  - `LAS::Property::gps_time` with type `double`
///  - `LAS::Property::R` with type `unsigned short`
///  - `LAS::Property::G` with type `unsigned short`
///  - `LAS::Property::B` with type `unsigned short`
///  - `LAS::Property::I` with type `unsigned short`
///
/// @sa `LAS::point_reader()`
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PropertyHandler handlers to recover properties.
///
/// @return true on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename ... PropertyHandler>
bool read_las_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef OutputIteratorValueType Enriched_point;

  LASreaderLAS lasreader;
  lasreader.open(stream);

  while (lasreader.read_point())
    {
      const LASpoint& laspoint = lasreader.point;
      Enriched_point new_point;

      LAS::internal::process_properties (laspoint, new_point, properties...);

      *(output ++) = new_point;
    }

  lasreader.close();
  
  return true;

}


/// \cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_las_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;

  return read_las_points_with_properties<OutputValueType>
    (stream, output, properties...);
}
/// \endcond

//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (position only) from a .las or .laz stream.
/// Potential additional properties are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value_type `CGAL::Point_3`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `CGAL::Point_3`.
///
/// @return `true` on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_las_points_with_properties (stream, output,
                                          LAS::point_reader(point_pmap));
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  // just deduce value_type of OutputIterator
  return read_las_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  return read_las_points
    <OutputIteratorValueType>(stream,
                              output,
                              make_identity_property_map(OutputIteratorValueType())
                              );
}

template < typename OutputIterator>
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  // just deduce value_type of OutputIterator
  return read_las_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output);
}
//-----------------------------------------------------------------------------------

/// @endcond


} //namespace CGAL

#endif // CGAL_READ_LAS_POINTS_H
