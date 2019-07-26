
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <CGAL/tags.h>

#include <iostream>


namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_, class Concurrency_tag_ = CGAL::Sequential_tag>
struct GarlandHeckbert_edge_collapse_visitor_base : Edge_collapse_visitor_base<TM_, Concurrency_tag_> {
  typedef Concurrency_tag_ Concurrency_tag;
  typedef TM_ TM;
  typedef typename Edge_collapse_visitor_base<TM>::Profile Profile;
  typedef typename Edge_collapse_visitor_base<TM>::vertex_descriptor vertex_descriptor;
  typedef typename Edge_collapse_visitor_base<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename Edge_collapse_visitor_base<TM>::Kernel Kernel;

  typedef typename internal::GarlandHeckbertCore<TM>::Matrix4x4 Matrix4x4;
  typedef typename internal::GarlandHeckbertCore<TM>::garland_heckbert_state_type garland_heckbert_state_type;
  typedef typename internal::GarlandHeckbertCore<TM>::FT                        FT;

  garland_heckbert_state_type& mCostMatrices;
  FT mDiscontinuityMultiplier;

  GarlandHeckbert_edge_collapse_visitor_base(
    garland_heckbert_state_type& aCostMatrices, FT aDiscontinuityMultiplier = 100.0)
    :
    mCostMatrices(aCostMatrices),
    mDiscontinuityMultiplier(aDiscontinuityMultiplier) {

  }


  void OnStarted(TM& aTM) {
    for(halfedge_descriptor& hd: halfedges(aTM)) {
      vertex_descriptor target_vertex = target(hd, aTM);
      if(mCostMatrices.find(target_vertex) == mCostMatrices.end()) {
        mCostMatrices.emplace(
          target_vertex,
          std::move(internal::GarlandHeckbertCore<TM>::fundamental_error_quadric(hd, aTM, mDiscontinuityMultiplier))
        );
      }
    }
  }

  void OnCollapsed(const Profile& aProfile,
    const vertex_descriptor& aVD) {


    Matrix4x4 combinedMatrix =
      std::move(internal::GarlandHeckbertCore<TM>::combine_matrices(
                          mCostMatrices.at(aProfile.v0()),
                          mCostMatrices.at(aProfile.v1())
                        ));

    mCostMatrices.erase(aProfile.v0());
    mCostMatrices.erase(aProfile.v1());

    mCostMatrices.emplace(aVD, std::move(combinedMatrix));

  }
};
}
}


#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H
