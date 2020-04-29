#ifndef CGAL_SHAPE_REGULARIZATION_TESTS_SAVER_H
#define CGAL_SHAPE_REGULARIZATION_TESTS_SAVER_H

// STL includes.
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// CGAL includes.
#include <CGAL/IO/io.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Shape_regularization {
namespace Tests {

  template<typename GeomTraits>
  struct Saver {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Polyline = std::vector<Point_3>;

    Saver() { 
      out.precision(20); 
    }
 
    inline std::string data() const {
      return out.str();
    }

    void export_polylines(
      const std::vector<Segment_2>& segments,
      const std::string file_path) {
      
      std::vector<Polyline> polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& s = segments[i].source();
        const auto& t = segments[i].target();
        
        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      export_polylines(polylines, file_path);
    }

    void export_group(
      const std::vector<Segment_2>& segments,
      const std::vector<std::size_t>& group,
      const std::string name) {
      
      std::vector<Segment_2> edges;
      for (const std::size_t seg_index : group)
        edges.push_back(segments[seg_index]);
      export_polylines(edges, "/Users/monet/Documents/gsoc/ggr/logs/" + name);
    }

    void export_closed_contour(
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

    void export_open_contour(
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
          "Error: Cannot save the file: " << path << std::endl; return;
      }
      
      file << data() << std::endl; file.close();
      std::cout << 
        "* data are saved in " << path << std::endl;
    }

    void export_polylines(
      const std::vector<Polyline>& polylines,
      const std::string file_path) {

      if (polylines.size() == 0)
        return;

      clear();
      for (std::size_t i = 0; i < polylines.size(); ++i) {
        const auto& polyline = polylines[i];

        out << polyline.size() << " ";
        for (std::size_t j = 0; j < polyline.size(); ++j)
          out << polyline[j] << " ";
        out << std::endl;
      }
      save(file_path + ".polylines");
    }
  };

} // namespace Tests
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_TESTS_SAVER_H
