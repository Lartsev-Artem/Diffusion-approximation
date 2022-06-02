#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include<vector>
#include<ctime>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/output_to_vtu.h>



// Domain

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT(Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
// Function
FT sphere_function(const Point& p)
{
    return CGAL::squared_distance(p, Point(0, 0, 0)) - pow(1./2, 2);
}
int maint()
{
    Mesh_domain domain =
        Mesh_domain::create_implicit_mesh_domain(sphere_function,
            K::Sphere_3(CGAL::ORIGIN, 2.));
    // Mesh criteria
    /*Mesh_criteria criteria(facet_angle = 10, facet_size = 0.02, facet_distance = 0.01,
        cell_radius_edge_ratio = 2, cell_size = 0.02);*/
    Mesh_criteria criteria(facet_angle = 30, facet_size = 0.1, facet_distance = 0.025,
        cell_radius_edge_ratio = 2, cell_size = 0.1);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    // Output
   // std::ofstream medit_file("C:\\Users\\Artem\\Desktop\\out");
    //c3t3.output_to_medit(medit_file);
    std::ofstream file("C:\\Users\\Artem\\Desktop\\TetraGrid\\Sphere.vtu");
   
    CGAL::output_to_vtu(file, c3t3, CGAL::IO::ASCII);
    return 0;
}

//    std::ofstream file("C:\\Users\\Artem\\Desktop\\out.vtu");
//    CGAL::output_to_vtu(file, c3t3);
//#include <CGAL/IO/output_to_vtu.h>
//
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Mesh_triangulation_3.h>
//#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
//#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/Labeled_mesh_domain_3.h>
//#include <CGAL/make_mesh_3.h>
//#include <CGAL/perturb_mesh_3.h>
//#include <CGAL/exude_mesh_3.h>
//// Domain
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef K::FT FT;
//typedef K::Point_3 Point;
//typedef FT(Function)(const Point&);
//typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
//// Triangulation
//typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
//typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
//// Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
//typedef Mesh_criteria::Facet_criteria Facet_criteria;
//typedef Mesh_criteria::Cell_criteria Cell_criteria;
//// Function
//FT ellipsoid_function(const Point& p)
//{
//    const FT x2 = p.x() * p.x();
//    const FT y2 = p.y() * p.y();
//    const FT z2 = p.z() * p.z();
//    return x2 + 2 * y2 + 4 * z2 - 1;
//}
//// To avoid verbose function and named parameters call
//using namespace CGAL::parameters;
//int main()
//{
//    // Domain (Warning: Sphere_3 constructor uses square radius !)
//    Mesh_domain domain =
//        Mesh_domain::create_implicit_mesh_domain(ellipsoid_function,
//            K::Sphere_3(CGAL::ORIGIN, 2.));
//    // Criteria
//    Facet_criteria facet_criteria(30, 0.08, 0.025); // angle, size, approximation
//    Cell_criteria cell_criteria(2, 0.1); // radius-edge ratio, size
//    Mesh_criteria criteria(facet_criteria, cell_criteria);
//    // Mesh generation (without optimization)
//    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
//    // Output
//    std::ofstream medit_file("C:\\Users\\Artem\\Desktop\\out_wo.vtu");
//    //c3t3.output_to_medit(medit_file)
//        CGAL::output_to_vtu(medit_file, c3t3,CGAL::IO::ASCII);
//    medit_file.close();
//    // Perturbation (5s, 12degree)
//    CGAL::perturb_mesh_3(c3t3, domain, time_limit = 5, sliver_bound = 12);
//    // Exudation
//    CGAL::exude_mesh_3(c3t3);
//    // Output
//    medit_file.open("C:\\Users\\Artem\\Desktop\\out_optimized.vtu");
//    CGAL::output_to_vtu(medit_file, c3t3);
//  //  c3t3.output_to_medit(medit_file);
//    return 0;
//}