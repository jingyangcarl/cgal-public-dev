// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_ANGLE_CONDITIONS_2_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_ANGLE_CONDITIONS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class Angle_conditions_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Angle_conditions_2() :
    m_margin_of_error(FT(1) / FT(4)) 
    { }

    FT reference(
      const Segment_data& seg_data, 
      const FT suffix) const {

      FT angle_deg = seg_data.orientation + suffix; 
      if (angle_deg < FT(0)) angle_deg += FT(180); 
      else if (angle_deg > FT(180)) angle_deg -= FT(180);
      return angle_deg;
    }

    int group_index(
      const FT in, const FT val_j, const int g_index) const {

      int g_j = -1;
      for (int k = -1; k <= 1; ++k) {  
        const FT tr = val_j - in + static_cast<FT>(k) * FT(180);
        if (CGAL::abs(tr) < m_margin_of_error) {  
          g_j = g_index; break;
        }
      }
      return g_j;
    }

    FT get_margin_of_error() const {
      CGAL_precondition(
        m_margin_of_error > FT(0));
      return m_margin_of_error;
    }

    void set_margin_of_error(const FT max_bound) {
      CGAL_precondition(max_bound > FT(0));
      m_margin_of_error = max_bound / FT(100);
    }

  private:
    FT m_margin_of_error;
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_ANGLE_CONDITIONS_2_H
