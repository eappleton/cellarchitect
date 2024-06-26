#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <iostream>
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  const char* fname = "examples/sphere1.off";
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();
  
  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }
  // Create domain
  Mesh_domain domain(polyhedron);
  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=40, facet_distance=20,
                         cell_radius_edge_ratio=2, cell_size=20);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, lloyd(time_limit=60), odt(time_limit=60), perturb(time_limit=30), exude(time_limit=30));
  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  printf("Done");
  return 0;
}