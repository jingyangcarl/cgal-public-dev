#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

using Kernel    = CGAL::Simple_cartesian<double>;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;

namespace CGAL {
namespace Shape_regularization {

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
class USER_quadratic_program  {
public:
  void reserve_d(const std::size_t) { }
  void reserve_c(const std::size_t) { }
  void reserve_a(const std::size_t) { }
  void reserve_b(const std::size_t) { }
  void reserve_l(const std::size_t) { }
  void reserve_u(const std::size_t) { }

  void  set_d(const std::size_t, const std::size_t, const FT) { }
  void  set_c(const std::size_t, const FT) { }
  void set_c0(const FT) { }
  void  set_a(const std::size_t, const std::size_t, const FT) { }
  void  set_b(const std::size_t, const FT) { }
  void  set_l(const std::size_t, const bool, const FT) { }
  void  set_u(const std::size_t, const bool, const FT) { }

  bool solve(std::vector<NT>& solution) {
    solution.clear();
    solution.resize(3, NT(0)); // 3 = 2 segments + 1 edge between them
    return true;
  }
};

} // namespace Shape_regularization
} // namespace CGAL

// Choose a type of a solver.
// #define OSQP_SOLVER
// #define CGAL_SOLVER
#define USER_SOLVER

#if defined(OSQP_SOLVER)
using Quadratic_program =
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>; // OSQP sparse solver
#endif
#if defined(CGAL_SOLVER)
using Quadratic_program =
  CGAL::Shape_regularization::CGAL_quadratic_program<FT>; // CGAL dense solver
#endif
#if defined(USER_SOLVER)
using Quadratic_program =
  CGAL::Shape_regularization::USER_quadratic_program<FT>; // USER custom solver
#endif

using NQ = CGAL::Shape_regularization::Custom_neighbor_query_2;
using RT = CGAL::Shape_regularization::Custom_regularization_2;
using QP = Quadratic_program;
using Regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, NQ, RT, QP>;

int main(int argc, char *argv[]) {

  NQ neighbor_query;
  RT angle_regularization;
  QP quadratic_program;

  Segments segments = {
    Segment_2(Point_2(-1,  0), Point_2(1, 0)),
    Segment_2(Point_2( 0, -1), Point_2(0, 1))
  };

  Regularizer regularizer(
    segments, neighbor_query, angle_regularization, quadratic_program, Kernel());
  regularizer.regularize();

  std::cout << "* regularized 2 segments" << std::endl;
}
