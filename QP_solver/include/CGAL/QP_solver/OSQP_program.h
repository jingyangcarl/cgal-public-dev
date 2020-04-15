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

#ifndef CGAL_OSQP_PROGRAM_H
#define CGAL_OSQP_PROGRAM_H

#include <CGAL/license/QP_solver.h>

// Warnings.
#include <CGAL/disable_warnings.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/functors.h>

// OSQP includes.
#include <osqp/osqp.h>

namespace CGAL {

  template<typename FT>
  class OSQP_program { 

  public:
    
    OSQP_program(
      CGAL::Comparison_result relation = CGAL::SMALLER,
      bool finite_lower = true,
      FT lower = 0,
      bool finite_upper = false,
      FT upper = 0) : 
    m_default_r(relation), 
    m_default_fl(finite_lower),
    m_default_l(lower), 
    m_default_fu(finite_upper), 
    m_default_u(upper), 
    m_is_valid(true) {
    
    CGAL_qpe_assertion(
      !finite_lower || !finite_upper || lower <= upper);
    }

  private:
    const CGAL::Comparison_result m_default_r;
    const bool m_default_fl, m_default_fu;
    const FT m_default_l, m_default_u;
    const bool m_is_valid;
  };

} // namespace CGAL

#endif // CGAL_OSQP_PROGRAM_H
