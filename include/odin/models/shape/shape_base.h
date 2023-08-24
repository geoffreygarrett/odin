#ifndef BASE_SHAPE_H
#define BASE_SHAPE_H

#include "include/tbb/parallel_for.h"
//#include "shape_concepts.h"
#include <Eigen/Core>
#include <array>
#include <odin/core/point.h>
#include <tbb/parallel_for.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <type_traits>

namespace odin::shape {

template<typename T, typename U, Point P, std::size_t dim>
class shape_base {
public:
    using derived_type      = T;
    using scalar_type       = U;
    using point_type        = P;// User-provided point container
    using point_series_type = typename point_series_type_trait<point_type>::type;
    //    using point_grid_type   = point_grid_type_trait<point_type>::type;

    // constructor
    explicit shape_base(point_type center = zero<point_type>())
        : m_center(center) {}

    // interface functions
    scalar_type volume() { return as_derived().volume_impl(); }

    scalar_type surface_area() { return as_derived().surface_area_impl(); }

    bool is_inside(const point_type &point) {
        return as_derived().is_inside_impl(point);
    }

    std::vector<bool> is_inside(const point_series_type &points) {
        size_t num_points = size_trait<point_series_type>::size(points);
        std::vector<bool> result(num_points);
        tbb::parallel_for(std::size_t(0), num_points, [&](std::size_t i) {
            result[i] = as_derived().is_inside_impl(
                    get_trait<point_series_type>::get(points, i));
        });
        return result;
    }

    point_type centroid() { return as_derived().centroid_impl(); }

    // getters and setters
    point_type get_center() const { return m_center; }

    void set_center(point_type center) { m_center = center; }

    // Interface function
    std::optional<point_type> ray_intersection(const point_type &origin,
                                               const point_type &direction) {
        // TODO: Consider https://www.geometrictools.com/index.html if
        //  we need more general ray intersection algorithms.
        // - https://en.wikipedia.org/wiki/CAP_theorem
        return as_derived().ray_intersection_impl(origin, direction);
    }

    //        point_series_type
    //        ray_intersection_series(const point_series_type &origins,
    //                                const point_series_type &directions) {
    //
    //            size_t num_points = origins.size();
    //
    //            thrust::device_vector<point_type> d_origins(
    //                    origins.data(), origins.data() + num_points);
    //            thrust::device_vector<point_type> d_directions(
    //                    directions.data(), directions.data() + num_points);
    //            thrust::device_vector<std::optional<point_type>> d_intersections(
    //                    num_points);
    //
    //            thrust::transform(d_origins.begin(), d_origins.end(),
    //                              d_directions.begin(), d_intersections.begin(),
    //                              [this] __device__(const point_type &origin,
    //                                                const point_type &direction) {
    //                                  return as_derived().ray_intersection_impl(
    //                                          origin, direction);
    //                              });
    //
    //            point_series_type intersections;
    //            intersections.reserve(num_points);
    //            thrust::copy_if(
    //                    thrust::device, d_intersections.begin(),
    //                    d_intersections.end(), std::back_inserter(intersections),
    //                    [] __device__(const std::optional<point_type> &opt) {
    //                        return opt.has_value();
    //                    });
    //
    //            return intersections;
    //        }


    point_series_type
    ray_intersection_series(const point_series_type &origins,
                            const point_series_type &directions) {

        size_t num_points = size_trait<point_series_type>::size(origins);

        // Preallocate space for intersections
        std::vector<std::optional<point_type>> intersections(num_points);

        // Compute ray-ellipsoid intersections in parallel
        tbb::parallel_for(
                tbb::blocked_range<size_t>(0, num_points),
                [&](const tbb::blocked_range<size_t> &r) {
                    for (size_t i = r.begin(); i != r.end(); ++i) {
                        auto origin
                                = get_trait<point_series_type>::get(origins, i);
                        auto direction = get_trait<point_series_type>::get(
                                directions, i);
                        auto intersection = as_derived().ray_intersection_impl(
                                origin, direction);
                        intersections[i] = intersection;
                    }
                });

        // Filter out invalid intersections
        auto is_invalid = [](const std::optional<point_type> &opt) {
            return !opt.has_value();
        };
        intersections.erase(std::remove_if(intersections.begin(),
                                           intersections.end(), is_invalid),
                            intersections.end());

        // Allocate space for valid intersections
        point_series_type valid_intersections;
        resize_trait<point_series_type>::resize(valid_intersections,
                                                intersections.size());

        // Populate valid_intersections with valid points
        size_t index = 0;
        for (const auto &intersection: intersections) {
            push_trait<point_series_type>::push(valid_intersections,
                                                *intersection, index);
            ++index;
        }

        return valid_intersections;
    }

    __device__ point_type invalid_point() {
        return point_type{std::numeric_limits<scalar_type>::quiet_NaN(),
                          std::numeric_limits<scalar_type>::quiet_NaN(),
                          std::numeric_limits<scalar_type>::quiet_NaN()};
    }

    static __device__ bool is_valid(const point_type &pt) {
        auto isnan = [](auto val) { return std::isnan(val); };
        return !std::any_of(pt.begin(), pt.end(), isnan);
    }

    thrust::device_vector<point_type> ray_intersection_series2(
            const thrust::device_vector<point_type> &origins,
            const thrust::device_vector<point_type> &directions) {
        assert(origins.size() == directions.size());
        size_t num_points = origins.size();

        // Create a vector to hold intersections
        thrust::device_vector<point_type> intersections(num_points);

        thrust::transform(thrust::device, origins.begin(), origins.end(),
                          directions.begin(), intersections.begin(),
                          [=] __device__(const point_type &origin,
                                         const point_type &direction) {
                              auto intersection
                                      = as_derived().ray_intersection_impl(
                                              origin, direction);
                              return intersection.has_value()
                                           ? intersection.value()
                                           : invalid_point();
                          });

        // Remove invalid intersections
        auto new_end = thrust::remove_if(
                thrust::device, intersections.begin(), intersections.end(),
                [](const point_type &pt) { return !is_valid(pt); });
        intersections.erase(new_end, intersections.end());

        return intersections;
    }

private:
    // private member functions
    derived_type &as_derived() { return static_cast<derived_type &>(*this); }

    // private data members
    P m_center;
};

}// namespace odin::shape

#endif//BASE_SHAPE_H
