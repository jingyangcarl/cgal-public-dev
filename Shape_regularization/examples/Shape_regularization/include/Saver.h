#ifndef CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H
#define CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// CGAL includes.
#include <CGAL/IO/io.h>

namespace CGAL {
namespace Shape_regularization {
namespace Examples {

  template<typename GeomTraits>
  struct Saver {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Polygon_3 = std::vector<Point_3>;

    inline std::string data() const {
      return out.str();
    }

    void save_segments_2(
      const std::vector<Segment_2>& segments, 
      const std::string path) {
      
      clear();
      for (const auto& segment : segments) {
        out << "v " << segment.source() << " " << FT(0) << std::endl;
        out << "v " << segment.target() << " " << FT(0) << std::endl;
        out << "v " << segment.target() << " " << FT(0) << std::endl;
      }
      for (std::size_t i = 0; i < segments.size() * 3; i += 3)
        out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << std::endl;
      save(path + ".obj");
    }

    void save_closed_contour_2(
      const std::vector<Point_2>& contour, 
      const std::string path) {

      clear();
    }

    void save_open_contour_2(
      const std::vector<Point_2>& contour, 
      const std::string path) {

      clear();
    }

    void save_polygons_3(
      const std::vector<Polygon_3>& polygons, 
      const std::string path) {

      clear();
    }

  private:
    std::stringstream out;

    void clear() {
      out.str(std::string());
    }

    void save(
      const std::string path) const {
      
      std::ofstream file(path.c_str(), std::ios_base::out);
      CGAL::set_ascii_mode(file);
      file.precision(20);
      if (!file) {
        std::cout << 
        "Error: cannot save the file: " << path << std::endl; return;
      }
      file << data() << std::endl;
      file.close();

      std::cout << 
        "* segments are saved in " 
      << path << std::endl;
    }
  };

} // namespace Examples
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H
