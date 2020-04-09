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
    using Polyline = std::vector<Point_3>;
    using Polylines = std::vector<Polyline>;

    Saver() { 
      out.precision(20); 
    }

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
      const std::string name) {

      if (contour.size() == 0)
        return;

      std::vector<Segment_2> segments;
      const std::size_t n = contour.size();
      segments.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;

        const auto& p = contour[i];
        const auto& q = contour[ip];
        segments.push_back(Segment_2(p, q));
      }
      export_polylines(segments, name);
    }

    void save_open_contour_2(
      const std::vector<Point_2>& contour, 
      const std::string name) {

      if (contour.size() == 0)
        return;

      std::vector<Segment_2> segments;
      const std::size_t n = contour.size();
      segments.reserve(n - 1);

      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& p = contour[i];
        const auto& q = contour[ip];
        segments.push_back(Segment_2(p, q));
      }
      export_polylines(segments, name);
    }

    template<typename Segment_wrapper_2>
    void export_polylines(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::string file_path) {

      std::vector<Segment_2> segments;
      segments.reserve(wraps.size());
      for (const auto& wrap : wraps)
        segments.push_back(wrap.segment);
      export_polylines(segments, file_path);
    }

    void export_polylines(
      const std::vector<Segment_2>& segments,
      const std::string file_path) {
      
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();
        
        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      export_polylines(polylines, file_path);
    }

    void export_polylines(
      const Polylines& polylines,
      const std::string file_path) {

      if (polylines.size() == 0)
        return;

      clear();
      for (std::size_t i = 0; i < polylines.size(); ++i) {
        const auto &polyline = polylines[i];

        out << polyline.size() << " ";
        for (std::size_t j = 0; j < polyline.size(); ++j)
          out << polyline[j] << " ";
        out << std::endl;
      }
      save(file_path + ".polylines");
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
