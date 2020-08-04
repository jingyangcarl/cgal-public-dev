#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

using Kernel    = CGAL::Simple_cartesian<double>;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;

namespace CGAL {

struct Custom_neighbor_query_2 {
  void operator()(
    const std::size_t query_index,
    std::vector<std::size_t>& neighbors) {
    neighbors.clear();
    if (query_index == 0) { neighbors.push_back(1); } // first  segment
    if (query_index == 1) { neighbors.push_back(0); } // second segment
  }
};

struct Custom_regularization_2 {
  FT bound(
    const std::size_t query_index) const {
    return FT(5); // max angle change
  }
  FT target(
    const std::size_t query_index_i,
    const std::size_t query_index_j) {
    return FT(0); // 0 angle change
  }
  void update(
    const std::vector<FT>& solution) {
    // skip update
  }
};

template<typename NT>
class USER_quadratic_program_traits  {
public:
  void reserve_P(const std::size_t) { }
  void reserve_q(const std::size_t) { }
  void reserve_A(const std::size_t) { }
  void reserve_l(const std::size_t) { }
  void reserve_u(const std::size_t) { }

  void  set_P(const std::size_t, const std::size_t, const FT) { }
  void  set_q(const std::size_t, const FT) { }
  void  set_r(const FT) { }
  void  set_A(const std::size_t, const std::size_t, const FT) { }
  void  set_l(const std::size_t, const FT) { }
  void  set_u(const std::size_t, const FT) { }

  template<typename OutputIterator>
  bool solve(OutputIterator solution) {

    // 3 = 2 segments + 1 edge between them
    for (std::size_t i = 0; i < 3; ++i)
      *(++solution) = NT(0);
    return true;
  }
};

} // namespace CGAL

// Choose a type of a solver.
// #define OSQP_SOLVER - be sure that OSQP is installed on your system!
#define USER_SOLVER

#if defined(OSQP_SOLVER)
using Quadratic_program =
  CGAL::OSQP_quadratic_program_traits<FT>; // OSQP sparse solver
#endif
#if defined(USER_SOLVER)
using Quadratic_program =
  CGAL::USER_quadratic_program_traits<FT>; // USER custom solver
#endif

using NQ = CGAL::Custom_neighbor_query_2;
using RT = CGAL::Custom_regularization_2;
using QP = Quadratic_program;
using Regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, NQ, RT, QP>;

int main(int argc, char *argv[]) {

  NQ neighbor_query;
  RT regularization_type;
  QP quadratic_program;

  Segments segments = {
    Segment_2(Point_2(-1,  0), Point_2(1, 0)),
    Segment_2(Point_2( 0, -1), Point_2(0, 1))
  };

  Regularizer regularizer(
    segments, neighbor_query, regularization_type, quadratic_program, Kernel());
  regularizer.regularize();

  std::cout << "* regularized 2 segments" << std::endl;
}
