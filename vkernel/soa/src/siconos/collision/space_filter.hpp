#pragma once

#include <cmath>
#include <concepts>
#include <numbers>
#include <print>
#include <stdexcept>

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/array_hash.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/diskfsegment_r.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/polytopedisk_r.hpp"
#include "siconos/collision/polytopefdisk_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
#include "siconos/collision/polytopemesh_r.hpp"
#include "siconos/collision/shape/chained_segment.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/mesh.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/collision/translated.hpp"
#include "siconos/simul/relation_sets.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace std {

// siconos::storage::mp::pair
template <typename First, typename Second>
struct hash<siconos::storage::mp::pair<First, Second>> {
  size_t operator()(const siconos::storage::mp::pair<First, Second>& p) const noexcept {
    size_t h1 = std::hash<First>{}(siconos::storage::mp::first(p));
    size_t h2 = std::hash<Second>{}(siconos::storage::mp::second(p));
    return h1 ^ (h2 + 0x9e3779b97f4a7c15 + (h1 << 6) + (h1 >> 2));
  }
};

// siconos::storage::mp::tuple<First, Second, Third>
template <typename First, typename Second, typename Third>
struct hash<siconos::storage::mp::tuple<First, Second, Third>> {
  size_t operator()(
      const siconos::storage::mp::tuple<First, Second, Third>& p) const noexcept {
    size_t h1 = std::hash<First>{}(siconos::storage::mp::tuple_first(p));
    size_t h2 = std::hash<Second>{}(siconos::storage::mp::tuple_second(p));
    size_t h3 = std::hash<Third>{}(siconos::storage::mp::tuple_third(p));
    size_t seed = 0;
    seed ^= h1 + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    seed ^= h2 + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    return seed;
  }
};
}  // namespace std

namespace siconos::collision {

template <typename Topology, typename Neighborhood>
struct space_filter : item {
  using items = gather<Topology, Neighborhood>;
  using topology = Topology;
  using dynamical_system = typename topology::fixed_dof_system;
  using finteraction = typename topology::finteraction;
  using dfinteraction = typename topology::dfinteraction;
  using nslaw = typename finteraction::nslaw;
  using neighborhood = Neighborhood;
  struct attributes {
    some::item_ref<topology> topology;
    some::item_ref<neighborhood> neighborhood;
    some::map<some::array<some::indice, some::indice_value<2>>, some::item_ref<nslaw>> nslaws;
    some::item_ref<diskdisk_r> diskdisk_r;
    some::map<some::array<some::indice, some::indice_value<2>>, some::item_ref<diskmesh_r>>
        diskmeshes;
    some::map<some::array<some::scalar, some::indice_value<4>>, some::item_ref<diskfsegment_r>>
        diskfsegments;
    // keyed by (segment x1 x2 y1 y2, manifold point index)
    some::map<some::array<some::scalar, some::indice_value<5>>,
              some::item_ref<polytopefsegment_r<Topology>>>
        polytopefsegments;
    some::map<some::array<some::scalar, some::indice_value<2>>,
              some::item_ref<collision::polytopefdisk_r<Topology>>>
        polytopefdisks;
    // keyed by (mesh, edge, manifold point index)
    some::map<some::array<some::indice, some::indice_value<3>>,
              some::item_ref<polytopemesh_r<Topology>>>
        polytopemeshes;
    some::map<some::array<some::scalar, some::indice_value<2>>, some::item_ref<diskfdisk_r>>
        diskfdisks;
    // the disk x polytope relation is stateless: one item, under key 0,
    // serves every mixed pair
    some::map<some::indice, some::item_ref<polytopedisk_r<Topology>>> polytopedisks;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaws() { return storage::attr<"nslaws">(*self()); }
    decltype(auto) insert_nonsmooth_law(auto nslaw, auto gid1, auto gid2) {
      auto& nslaws = storage::attr<"nslaws">(*self());
      nslaws[{static_cast<std::size_t>(gid1), static_cast<std::size_t>(gid2)}] = nslaw;
      nslaws[{static_cast<std::size_t>(gid2), static_cast<std::size_t>(gid1)}] = nslaw;
      return *self();
    }
    // the law of the collision groups of two bodies (or body and static shape)
    decltype(auto) nonsmooth_law(auto& body1, auto& body2) {
      auto gid1 = static_cast<std::size_t>(body1.collision_group());
      auto gid2 = static_cast<std::size_t>(body2.collision_group());
      auto& nslaws = storage::attr<"nslaws">(*self());
      auto search = nslaws.find({gid1, gid2});
      if (search == nslaws.end()) {
        throw std::runtime_error(
            "space_filter: no nonsmooth law registered for collision groups (" +
            std::to_string(gid1) + ", " + std::to_string(gid2) + ")");
      }
      return storage::make_handle(self()->data(), search->second);
    }
    decltype(auto) topology() {
      return storage::make_ref_handle(self()->data(), storage::attr<"topology">(*self()));
    }

    decltype(auto) neighborhood() {
      return storage::make_ref_handle(self()->data(), storage::attr<"neighborhood">(*self()));
    }

    decltype(auto) diskdisk_r() {
      return storage::make_ref_handle(self()->data(), storage::attr<"diskdisk_r">(*self()));
    }

    decltype(auto) diskmeshes() { return storage::attr<"diskmeshes">(*self()); }

    decltype(auto) diskfsegments() { return storage::attr<"diskfsegments">(*self()); }

    decltype(auto) polytopefsegments() { return storage::attr<"polytopefsegments">(*self()); }

    decltype(auto) polytopefdisks() { return storage::attr<"polytopefdisks">(*self()); }

    decltype(auto) polytopemeshes() { return storage::attr<"polytopemeshes">(*self()); }

    decltype(auto) diskfdisks() { return storage::attr<"diskfdisks">(*self()); }
    decltype(auto) polytopedisks() { return storage::attr<"polytopedisks">(*self()); }

    // void create_relations()
    // {
    //   using env_t = decltype(self()->env());
    //   using indice = typename env_t::indice;

    //   using points_t = typename decltype(self()->neighborhood())::points_t;

    //   auto& data = self()->data();

    //   mp::for_each(points_t{}, []<typename Point>(Point) {
    //     using item_t = typename Point::item_t;
    //     auto& item_handles = storage::handles<item_t>(data);

    //     if constexpr (std::derived_from<item_t, shape::segment>) {
    //       for (auto segment_handle : item_handles) {
    //         auto rel = storage::add<diskfsegment_r>(data);
    //         rel.shape() = segment_handle;
    //       }
    //     }
    //     else if constexpr (std::derived_from<item_t,
    //                                          translated<shape::disk>>) {
    //       for (auto fdisk_handle : item_handles) {
    //         auto rel = storage::add<diskfdisk_r>(data);
    //         rel.shape() = fdisk_handle;
    //       }
    //     }
    //     else if constexpr (std::derived_from<item_t,
    //                                          shape::chained_segment>) {
    //       for (auto mesh_handle : item_handles) {
    //         indice nbsegments = std::size(mesh_handle.points()) / 2;
    //         for (indice i = 0; i < nbsegments; ++i) {
    //           auto rel = storage::add<diskmesh_r>(data);
    //           rel.shape() = mesh_handle;
    //           rel.contact_index() = i;
    //         }
    //       }
    //     }
    //   });
    // }

    // void __init__() { create_relations(); }

    void remove_static_item(auto step, auto item_handle) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using item_t = typename std::decay_t<decltype(item_handle)>::type;
      using point_t = collision::point<empty_item, item_t>;
      using points_t = typename neighborhood::points_t;

      auto& data = self()->data();

      auto& points_flags = storage::attr_values<point_t, "flag">(data, step);
      auto& points_items = storage::attr_values<point_t, "item">(data, step);
      auto& points_coords = storage::attr_values<point_t, "coord">(data, step);

      auto ps_indx = mp::index_of<point_t>(mp::std_tuple(points_t{}));

      // mp::index_if(
      //   points_t{},
      //   mp::equal.to(collision::point<item_t>));

      // first remove item
      storage::remove(data, item_handle);

      // find all associated points
      for (auto [flag] : view::zip(points_flags)) {
        flag = false;
      }

      for (auto [flag, pitem] : view::zip(points_flags, points_items)) {
        if (item_handle.index().value() == pitem.value()) {
          flag = true;
        };
      }
      // XX check with one loop : flag = (item_handle.get() == pitem.get())
      // remove points
      auto ff = std::ranges::find(points_flags, true);

      while (ff != points_flags.end()) {
        auto ff_index = ff - points_flags.begin();
        auto point = storage::make_handle(data, storage::index<point_t, indice>(ff_index));
        storage::remove(data, point);

        auto remaining_points_flags = std::ranges::subrange(ff, points_flags.end());
        auto rff = std::ranges::find(remaining_points_flags, true);

        ff = ff + (rff - remaining_points_flags.begin());
      };

      neighborhood().instance()->resize_point_set(ps_indx, points_coords.front().data(),
                                                  points_coords.size());
    }

    void insert_diskfsegment_r(auto hdl) {
      storage::attr<"diskfsegments">(*self())[{hdl.segment().x1(), hdl.segment().x2(),
                                               hdl.segment().y1(), hdl.segment().y2()}] = hdl;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void insert_polytopefsegment_r(auto hdl) {
      using scalar = std::decay_t<decltype(hdl.segment().x1())>;
      storage::attr<"polytopefsegments">(
          *self())[{hdl.segment().x1(), hdl.segment().x2(), hdl.segment().y1(),
                    hdl.segment().y2(), static_cast<scalar>(hdl.point_index())}] = hdl;
    }

    void insert_polytopemesh_r(auto hdl) {
      storage::attr<"polytopemeshes">(
          *self())[{hdl.mesh().index().value(), hdl.contact_index(), hdl.point_index()}] = hdl;
    }

    void insert_diskfdisk_r(auto tds) {
      auto htds = tds.translated_disk_shape();
      storage::attr<"diskfdisks">(*self())[{htds.translation()[0], htds.translation()[1]}] =
          tds;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void make_points() {
      auto& data = self()->data();
      using env = decltype(self()->env());
      using indice = typename env::indice;

      using points_t = typename neighborhood::points_t;

      mp::for_each(points_t{}, [&]<typename Point>(Point) {
        using item_t = typename Point::item_t;

        if constexpr (std::derived_from<item_t, model::lagrangian_ds> ||
                      std::derived_from<item_t, model::elastic_lagrangian_ds>) {
          auto all_ds = storage::handles<item_t>(data);
          for (auto ds : all_ds) {
            //            print("add disk point : {}\n", ds.get());
            if constexpr (collision::always_shape<collision::shape::mesh, decltype(ds)>()) {
              auto shape_handle = collision::shape_of<collision::shape::mesh>(ds);

              auto nbsegments = std::size(shape_handle.segments().nodes()) - 1;
              // multiple points are associated to the system
              for (std::size_t index = 0; index < nbsegments; ++index) {
                for (auto [i, point_coord] :
                     shape_handle.segments().points_coords(index) | view::enumerate) {
                  auto new_point = storage::add<Point>(data);
                  new_point.item() = ds;
                  // the mesh of this system: a contact on this point is on
                  // one of its edges
                  new_point.shape() = shape_handle;
                  new_point.coord() = algebra::cast(mp::type_c<float>, point_coord);
                  new_point.seg_index() = index;
                  new_point.point_index() = i;
                }
              }
            } else {
              // With one point at its center, a disk pair is seen by
              // the neighborhood search iff r1 + r2 <= search radius,
              // so a single point is only enough for 2 r <= search
              // radius. A bigger disk gets points on its circle, spaced
              // by half the search radius: any contact point is then
              // within a quarter of the search radius of one of them.
              // Polytopes keep one point at their center.
              auto npoints = std::size_t{1};
              auto radius = 0.;
              if constexpr (collision::may_have_shape<collision::shape::disk,
                                                      decltype(ds)>()) {
                if (collision::has_shape<collision::shape::disk>(ds)) {
                  auto search_radius = neighborhood().instance()->radius();
                  radius = collision::shape_of<collision::shape::disk>(ds).radius();
                  if (2 * radius > search_radius) {
                    npoints = static_cast<std::size_t>(
                        std::ceil(2 * std::numbers::pi * radius / (search_radius / 2)));
                  }
                }
              }
              for (std::size_t i = 0; i < npoints; ++i) {
                auto new_point = storage::add<Point>(data);
                new_point.item() = ds;
                if (npoints == 1) {
                  new_point.offset() = {0.f, 0.f};
                } else {
                  auto t = 2 * std::numbers::pi * i / npoints;
                  new_point.offset() = {static_cast<float>(radius * std::cos(t)),
                                        static_cast<float>(radius * std::sin(t))};
                }
                new_point.update(0);
              }
            }
          }
        }
        // fixed segment
        else if constexpr (std::derived_from<item_t, collision::shape::segment>) {
          auto all_segments = storage::handles<item_t>(data);
          for (auto segment : all_segments) {
            for (auto [i, point_coord] : segment.points_coords() | view::enumerate) {
              auto new_point = storage::add<Point>(data);
              new_point.item() = segment;
              new_point.coord() = algebra::cast(mp::type_c<float>, point_coord);
              new_point.point_index() = i;
            }
          }
        }
        // fixed disk
        else if constexpr (std::derived_from<item_t,
                                             collision::translated<collision::shape::disk>>) {
          auto all_fdisks = storage::handles<item_t>(data);
          for (auto fdisk : all_fdisks) {
            for (auto [i, point_coord] : fdisk.points_coords() | view::enumerate) {
              auto new_point = storage::add<Point>(data);
              new_point.item() = fdisk;
              new_point.coord() = algebra::cast(mp::type_c<float>, point_coord);
              new_point.point_index() = (indice)i;
            }
          }
        };
      });
    }

    template <typename Index1, typename Index2>
    decltype(auto) make_ipair(Index1 ids1, Index2 ids2) {
      auto i1 = ids1.value();
      auto i2 = ids2.value();

      if constexpr (std::is_same_v<Index1, Index2>)
      // same type of indices, relation is the same with permutation
      {
        if (i1 < i2) {
          return mp::make_pair(i1, i2);
        } else {
          return mp::make_pair(i2, i1);
        }
      } else {
        // not same type of indices: no permutation!
        return mp::make_pair(i1, i2);
      }
    }

    /* key of a ds-ds interaction: the ds pair (ordered as make_ipair) and
       the manifold point the interaction stands for. Polytope relations
       carry a point_index (a polytope pair has one interaction per point of
       its contact manifold); the others stand for the whole pair: 0. */
    template <typename Index1, typename Index2>
    auto make_ikey(Index1 ids1, Index2 ids2, std::size_t pi) {
      auto i1 = ids1.value();
      auto i2 = ids2.value();
      if constexpr (std::is_same_v<Index1, Index2>) {
        if (i2 < i1) std::swap(i1, i2);
      }
      return mp::make_tuple(i1, i2, static_cast<decltype(i1)>(pi));
    }

    std::size_t point_of(auto inter) {
      return siconos::variant::visit(self()->data(), inter.relation(),
                                     [](auto rrel) -> std::size_t {
                                       if constexpr (requires { rrel.point_index(); }) {
                                         return rrel.point_index();
                                       } else {
                                         return 0;
                                       }
                                     });
    }

    void build_proximity_maps(auto& ds_ds_prox, auto& ds_segment_prox, auto& ds_fdisk_prox,
                              const auto& ds1s, const auto& ds2s, const auto& interactions) {
      auto& data = self()->data();

      // build (ds ds -> inter) & (ds segment -> inter) maps
      for (auto [ds1, ds2, inter] : view::zip(ds1s, ds2s, interactions)) {
        if (ds1 != ds2) {
          ds_ds_prox[make_ikey(ds1, ds2, point_of(inter))] = inter;
        } else {
          siconos::variant::visit(
              data, inter.relation(),
              mp::overload(
                  // https://stackoverflow.com/questions/46114214/lambda-implicit-capture-fails-with-variable-declared-from-structured-binding
                  // capture by value ok with handles
                  [&, ds1 = ds1, inter = inter]<match::handle<diskfsegment_r> DiskSegmentR>(
                      DiskSegmentR rel) {
                    auto segment = rel.segment();
                    auto coefs =
                        std::array{segment.x1(), segment.x2(), segment.y1(), segment.y2(),
                                   std::decay_t<decltype(segment.x1())>{0}};
                    ds_segment_prox[mp::make_pair(ds1.value(), coefs)] = inter;
                  },
                  [&, ds1 = ds1,
                   inter =
                       inter]<match::handle<polytopefsegment_r<Topology>> PolytopeFsegmentR>(
                      PolytopeFsegmentR rel) {
                    auto segment = rel.segment();
                    auto coefs = std::array{
                        segment.x1(), segment.x2(), segment.y1(), segment.y2(),
                        static_cast<std::decay_t<decltype(segment.x1())>>(rel.point_index())};
                    ds_segment_prox[mp::make_pair(ds1.value(), coefs)] = inter;
                  },
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<diskfdisk_r> DiskFdiskR>(DiskFdiskR rel) {
                    auto fdisk = rel.translated_disk_shape();
                    auto& translat = fdisk.translation();
                    auto coefs = std::array{translat[0], translat[1]};
                    ds_fdisk_prox[mp::make_pair(ds1.value(), coefs)] = inter;
                  },
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<polytopefdisk_r<Topology>> PolytopeFdiskR>(
                      PolytopeFdiskR rel) {
                    auto& translat = rel.translated_disk_shape().translation();
                    auto coefs = std::array{translat[0], translat[1]};
                    ds_fdisk_prox[mp::make_pair(ds1.value(), coefs)] = inter;
                  },
                  []<bool flag = false>(auto) {
                    // hardening: this used to be a silent no-op in release
                    // (assert compiled out); a self interaction must bind a
                    // fixed-object relation
                    throw std::runtime_error(
                        "space_filter::build_proximity_maps: unhandled relation "
                        "for a self interaction (ds1 == ds2); only "
                        "diskfsegment_r, diskfdisk_r, polytopefsegment_r and polytopefdisk_r "
                        "are supported");
                  }));
        }
      }
    }

    /* third slot of a ds-mesh proximity key: the mesh edge and the manifold
       point of the interaction on it (a polytope has one interaction per
       point, two in 2D; a disk has one, point 0) */
    static auto edge_point_slot(auto edge, std::size_t pi) {
      return static_cast<decltype(edge)>(2 * edge + pi);
    }

    void build_dproximity_maps(auto& ds_dds_prox, const auto& ds1s, const auto& ds2s,
                               const auto& dinteractions) {
      auto& data = self()->data();

      // build (ds ds -> inter) & (ds segment -> inter) maps
      for (auto [ds1, ds2, inter] : view::zip(ds1s, ds2s, dinteractions)) {
        auto contact_index = variant::visit(data, inter.relation(),
                                            [](auto& rrel) { return rrel.contact_index(); });
        ds_dds_prox[mp::make_tuple(ds1.value(), ds2.value(),
                                   edge_point_slot(contact_index, point_of(inter)))] = inter;
      }
    }

    void dsds_activation(auto& ds1, auto& ds2, auto& fmap, auto hp1, auto hp2) {
      // two points of the same big disk are neighbors of each other
      if constexpr (std::is_same_v<std::decay_t<decltype(ds1)>, std::decay_t<decltype(ds2)>>) {
        if (ds1.index() == ds2.index()) return;
      }
      auto& data = self()->data();
      auto topo = self()->topology();
      auto nslaw = nonsmooth_law(ds1, ds2);
      auto diskdisk_r = self()->diskdisk_r();

      const bool polytope1 = collision::has_shape<collision::shape::polytope>(ds1);
      const bool polytope2 = collision::has_shape<collision::shape::polytope>(ds2);

      // one edge per disk pair; one per manifold point (two in 2D) for a
      // polytope pair, which the cxd stage rebinds to polytopepolytope_r
      const std::size_t npoints = polytope1 && polytope2 ? 2 : 1;
      for (std::size_t pi = 0; pi < npoints; ++pi) {
        auto find_inter = fmap.find(make_ikey(ds1.index(), ds2.index(), pi));
        if (find_inter != fmap.end()) {
          // keep this edge
          auto inter = storage::make_handle(data, std::get<1>(*find_inter));
          storage::prop<"activation", Topology>(inter) = true;
          continue;
        }
        if (polytope1 != polytope2) {
          // disk x polytope: linked disk first, as polytopedisk_r expects
          using ds_t = std::decay_t<decltype(ds1)>;
          using pd_index_t = storage::index<collision::polytopedisk_r<Topology>, std::size_t>;
          using inter_t = decltype(topo.link(ds1, ds2));
          if constexpr (collision::may_have_shape<collision::shape::polytope, ds_t>() &&
                        collision::may_have_shape<collision::shape::disk, ds_t>()) {
            if constexpr (requires(
                              std::decay_t<decltype(std::declval<inter_t&>().relation())> v,
                              pd_index_t i) { v = i; }) {
              auto inter = polytope1 ? topo.link(ds2, ds1) : topo.link(ds1, ds2);
              inter.nslaw() = nslaw;
              storage::prop<"activation", Topology>(inter) = true;
              auto& rmap = self()->polytopedisks();
              if (auto search = rmap.find(0); search != rmap.end()) {
                simul::bind_relation<Topology>(data, inter, search->second);
              } else {
                auto rel = storage::add<collision::polytopedisk_r<Topology>>(data);
                simul::bind_relation<Topology>(data, inter, rel);
                rmap[0] = rel;
              }
              fmap[make_ikey(ds1.index(), ds2.index(), pi)] = inter;
            } else {
              throw std::runtime_error(
                  "space_filter: a disk and a polytope are neighbors, but the interaction "
                  "relations of this world do not include polytopedisk_r");
            }
          }
          continue;
        }

        // the generic binder serves disk x disk, and polytope x polytope
        // until the cxd stage rebinds it
        // create the edge
        auto inter = topo.link(ds1, ds2);
        inter.nslaw() = nslaw;  // one nslaw for the moment

        storage::prop<"activation", Topology>(inter) = true;

        simul::bind_relation<Topology>(data, inter, diskdisk_r);
        fmap[make_ikey(ds1.index(), ds2.index(), pi)] = inter;
      }
    }

    void dsdds_activation(auto& ds1, auto& ds2, auto& dmap, auto hp1, auto hp2) {
      auto& data = self()->data();
      auto topo = self()->topology();
      auto nslaw = nonsmooth_law(ds1, ds2);
      auto& diskmeshes = self()->diskmeshes();
      auto& polytopemeshes = self()->polytopemeshes();

      auto mesh = hp2.shape();
      auto contact_index = hp2.seg_index();

      // the shape kind decides the relation: compile time for a static shape
      // declaration, a set membership for per-kind shapes. Decided before
      // the edge exists, so a refused pair leaves no interaction behind.
      using ds1_t = std::decay_t<decltype(ds1)>;
      const bool polytope1 = collision::has_shape<collision::shape::polytope>(ds1);
      const bool disk1 = collision::has_shape<collision::shape::disk>(ds1);
      if (!polytope1 && !disk1) {
        throw std::runtime_error(
            "diskmesh_r/polytopemesh_r requires a disk- or polytope-shaped ds; "
            "no narrowphase relation exists for this shape pair");
      }

      // one edge per (ds, mesh edge) for a disk, one per manifold point (two
      // in 2D) for a polytope
      const std::size_t npoints = polytope1 ? 2 : 1;
      for (std::size_t pi = 0; pi < npoints; ++pi) {
        const auto key = mp::make_tuple(ds1.index().value(), ds2.index().value(),
                                        edge_point_slot(contact_index, pi));
        if (auto find_inter = dmap.find(key); find_inter != dmap.end()) {
          // keep this edge
          auto inter = storage::make_handle(data, std::get<1>(*find_inter));
          storage::prop<"activation", Topology>(inter) = true;
          continue;
        }

        // create the edge
        auto inter = topo.link(ds1, ds2);
        inter.nslaw() = nslaw;  // one nslaw for the moment
        storage::prop<"activation", Topology>(inter) = true;

        if constexpr (collision::may_have_shape<collision::shape::polytope, ds1_t>()) {
          if (polytope1) {
            // polytope x mesh edge, manifold point pi
            using indice = std::decay_t<decltype(mesh.index().value())>;
            const std::array<indice, 3> rkey{mesh.index().value(),
                                             static_cast<indice>(contact_index),
                                             static_cast<indice>(pi)};
            if (auto search = polytopemeshes.find(rkey); search != polytopemeshes.end()) {
              simul::bind_relation<Topology>(data, inter, search->second);
            } else {
              auto rel = storage::add<collision::polytopemesh_r<Topology>>(data);
              rel.mesh() = mesh;
              rel.contact_index() = contact_index;
              rel.point_index() = pi;
              simul::bind_relation<Topology>(data, inter, rel);
              polytopemeshes[rkey] = rel;
            }
          }
        }
        if constexpr (collision::may_have_shape<collision::shape::disk, ds1_t>()) {
          if (disk1) {
            // disk x mesh edge
            if (auto search = diskmeshes.find({mesh.index().value(), contact_index});
                search != diskmeshes.end()) {
              simul::bind_relation<Topology>(data, inter, search->second);
            } else {
              auto rel = storage::add<collision::diskmesh_r>(data);
              rel.mesh() = mesh;
              rel.contact_index() = contact_index;
              simul::bind_relation<Topology>(data, inter, rel);
              diskmeshes[{mesh.index().value(), contact_index}] = rel;
            }
          }
        }
        dmap[key] = inter;
      }
    }

    void diskfsegment_activation(auto& body1, auto& body2, auto& intermap) {
      auto& data = self()->data();
      auto nslaw = nonsmooth_law(body1, body2);
      auto topo = self()->topology();
      auto& relmap = self()->diskfsegments();

      // diskfsegment_r reads shape.radius() on the body, which is correct only
      // for a disk: a polytope exposes radius() as its circumradius.
      if (!collision::has_shape<collision::shape::disk>(body1)) {
        throw std::runtime_error(
            "diskfsegment_r requires a disk-shaped ds; no narrowphase relation "
            "exists for this shape pair");
      }

      auto segment = body2;
      auto x1 = segment.x1();
      auto x2 = segment.x2();
      auto y1 = segment.y1();
      auto y2 = segment.y2();

      auto& ds1 = body1;

      auto find_inter = intermap.find(mp::make_pair(
          ds1.index().value(), std::array{x1, x2, y1, y2, std::decay_t<decltype(x1)>{0}}));

      if (find_inter != intermap.end()) {
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        // keep this interaction
        storage::prop<"activation", Topology>(inter) = true;

        // print("interaction FOUND for
        // {},{},{}\n", a, b,
        //       c);
      } else {
        // print("interaction NOT FOUND for
        // {},{},{}\n", a,
        //       b, c);
        // create the edge
        auto inter = topo.link(body1);
        inter.nslaw() = nslaw;  // one nslaw for the moment

        if (auto search = relmap.find({x1, x2, y1, y2}); search != relmap.end()) {
          simul::bind_relation<Topology>(data, inter, search->second);
        } else {
          auto dl = storage::add<diskfsegment_r>(data);

          // set segment pointer toward body2
          dl.segment() = segment;

          simul::bind_relation<Topology>(data, inter, dl);
          relmap[{x1, x2, y1, y2}] = dl;
        }
        storage::prop<"activation", Topology>(inter) = true;
        intermap[mp::make_pair(ds1.index().value(),
                               std::array{x1, x2, y1, y2, std::decay_t<decltype(x1)>{0}})] =
            inter;
      }
    }

    void polytopefsegment_activation(auto& body1, auto& body2, auto& intermap) {
      auto& data = self()->data();
      auto nslaw = nonsmooth_law(body1, body2);
      auto topo = self()->topology();
      auto& relmap = self()->polytopefsegments();

      if (!collision::has_shape<collision::shape::polytope>(body1)) {
        throw std::runtime_error(
            "polytopefsegment_r requires a polytope-shaped ds; no narrowphase "
            "relation exists for this shape pair");
      }

      auto segment = body2;
      auto x1 = segment.x1();
      auto x2 = segment.x2();
      auto y1 = segment.y1();
      auto y2 = segment.y2();
      using scalar = std::decay_t<decltype(x1)>;

      auto& ds1 = body1;

      // one interaction per point of the contact manifold (at most two in
      // 2D); a point the manifold does not have keeps an inactive gap
      using poly_fs_t = collision::polytopefsegment_r<Topology>;
      for (std::size_t pi = 0; pi < 2; ++pi) {
        const std::array key{x1, x2, y1, y2, static_cast<scalar>(pi)};
        auto find_inter = intermap.find(mp::make_pair(ds1.index().value(), key));

        if (find_inter != intermap.end()) {
          auto inter = storage::make_handle(data, std::get<1>(*find_inter));
          storage::prop<"activation", Topology>(inter) = true;
        } else {
          auto inter = topo.link(body1);
          inter.nslaw() = nslaw;

          if (auto search = relmap.find(key); search != relmap.end()) {
            simul::bind_relation<Topology>(data, inter, search->second);
          } else {
            auto dl = storage::add<poly_fs_t>(data);
            dl.segment() = segment;
            dl.point_index() = pi;
            simul::bind_relation<Topology>(data, inter, dl);
            relmap[key] = dl;
          }
          storage::prop<"activation", Topology>(inter) = true;
          intermap[mp::make_pair(ds1.index().value(), key)] = inter;
        }
      }
    }

    void polytopefdisk_activation(auto& body1, auto& body2, auto& intermap) {
      auto& data = self()->data();
      auto nslaw = nonsmooth_law(body1, body2);
      auto topo = self()->topology();
      auto& relmap = self()->polytopefdisks();

      if (!collision::has_shape<collision::shape::polytope>(body1)) {
        throw std::runtime_error(
            "polytopefdisk_r requires a polytope-shaped ds; no narrowphase "
            "relation exists for this shape pair");
      }

      auto fdisk = body2;
      auto& translat = fdisk.translation();
      auto coefs = std::array{translat[0], translat[1]};

      auto& ds1 = body1;

      auto find_inter = intermap.find(mp::make_pair(ds1.index().value(), coefs));

      if (find_inter != intermap.end()) {
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        storage::prop<"activation", Topology>(inter) = true;
      } else {
        auto inter = topo.link(body1);
        inter.nslaw() = nslaw;

        using poly_fd_t = collision::polytopefdisk_r<Topology>;
        if (auto search = relmap.find(coefs); search != relmap.end()) {
          simul::bind_relation<Topology>(data, inter, search->second);
        } else {
          auto dl = storage::add<poly_fd_t>(data);
          dl.translated_disk_shape() = fdisk;
          simul::bind_relation<Topology>(data, inter, dl);
          relmap[coefs] = dl;
        }
        storage::prop<"activation", Topology>(inter) = true;
        intermap[mp::make_pair(ds1.index().value(), coefs)] = inter;
      }
    }

    void diskfdisk_activation(auto& body1, auto& body2, auto& ds_fdisk_prox) {
      auto& data = self()->data();
      auto nslaw = nonsmooth_law(body1, body2);
      auto topo = self()->topology();
      auto& diskfdisks = self()->diskfdisks();

      auto fdisk = body2;
      auto& translat = fdisk.translation();
      auto coefs = std::array{translat[0], translat[1]};

      // diskfdisk_r reads shape.radius() on the body, which is correct only for
      // a disk: a polytope exposes radius() as its circumradius.
      if (!collision::has_shape<collision::shape::disk>(body1)) {
        throw std::runtime_error(
            "diskfdisk_r requires a disk-shaped ds; no narrowphase relation "
            "exists for this shape pair");
      }

      auto& ds1 = body1;

      auto find_inter = ds_fdisk_prox.find(mp::make_pair(ds1.index().value(), coefs));

      if (find_inter != ds_fdisk_prox.end()) {
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        storage::prop<"activation", Topology>(inter) = true;
      } else {
        auto inter = topo.link(body1);
        inter.nslaw() = nslaw;

        if (auto search = diskfdisks.find(coefs); search != diskfdisks.end()) {
          simul::bind_relation<Topology>(data, inter, search->second);
        } else {
          auto dfd = storage::add<diskfdisk_r>(data);

          // set disk pointer toward body2
          dfd.translated_disk_shape() = fdisk;

          simul::bind_relation<Topology>(data, inter, dfd);
          diskfdisks[coefs] = dfd;
        }
        storage::prop<"activation", Topology>(inter) = true;
        ds_fdisk_prox[mp::make_pair(ds1.index().value(), coefs)] = inter;
      }
    }

    template <typename Interaction>
    void remove_interactions(auto& ds_ds_prox, auto& ds_segment_prox, auto& ds_fdisk_prox) {
      using env = decltype(self()->env());
      using indice = typename env::indice;

      auto& data = self()->data();

      auto& activations = storage::prop_values<Interaction, "activation", Topology>(data, 0);

      auto interactions = storage::handles<Interaction>(data, 0);

      //      print("BEFORE REMOVAL: size of indexset0: {}\n",
      //      activations.size());
      for (auto [activation, inter] : view::zip(activations, interactions)) {
        if (!activation) {
          //          print("START REMOVE interaction {}\n", inter.get());

          if (storage::prop<"ds1", Topology>(inter) != storage::prop<"ds2", Topology>(inter)) {
            auto finter = ds_ds_prox.find(make_ikey(storage::prop<"ds1", Topology>(inter),
                                                    storage::prop<"ds2", Topology>(inter),
                                                    point_of(inter)));

            // hardening: duplicate ds pairs (the proximity map keeps only
            // the last interaction per pair) or a map rebuilt out of sync
            // would leave this find empty; the assert is compiled out in a
            // release build, so this used to dereference an end iterator
            if (finter == ds_ds_prox.end()) {
              throw std::runtime_error(
                  "space_filter::remove_interactions: no proximity map entry for "
                  "this ds-ds interaction (stale filter state or duplicate ds "
                  "pair)");
            }
            if (inter.index().value() != std::get<1>(*finter).value()) {
              throw std::runtime_error(
                  "space_filter::remove_interactions: two interactions share the "
                  "same ds pair; this one is stale");
            }

            ds_ds_prox.erase(finter);
            //            print("  REMOVE ds ds interaction between {}
            //            {}\n",
            //                  storage::prop<"ds1", Topology>(inter).get(),
            //                  storage::prop<"ds2", Topology>(inter).get());
          } else {
            siconos::variant::visit(
                data, inter.relation(),
                mp::overload(
                    [&]<match::handle<diskfsegment_r> DiskSegmentR>(DiskSegmentR rel) {
                      auto segment = rel.segment();
                      auto coefs =
                          std::array{segment.x1(), segment.x2(), segment.y1(), segment.y2(),
                                     std::decay_t<decltype(segment.x1())>{0}};
                      auto finter = ds_segment_prox.find(
                          mp::make_pair(storage::prop<"ds1", Topology>(inter).value(), coefs));
                      if (finter == ds_segment_prox.end()) {
                        throw std::runtime_error(
                            "space_filter::remove_interactions: no proximity map "
                            "entry for this ds-fsegment interaction");
                      }
                      ds_segment_prox.erase(finter);
                    },
                    [&]<match::handle<polytopefsegment_r<Topology>> PolytopeFsegmentR>(
                        PolytopeFsegmentR rel) {
                      auto segment = rel.segment();
                      auto coefs =
                          std::array{segment.x1(), segment.x2(), segment.y1(), segment.y2(),
                                     static_cast<std::decay_t<decltype(segment.x1())>>(
                                         rel.point_index())};
                      auto finter = ds_segment_prox.find(
                          mp::make_pair(storage::prop<"ds1", Topology>(inter).value(), coefs));
                      if (finter == ds_segment_prox.end()) {
                        throw std::runtime_error(
                            "space_filter::remove_interactions: no proximity map "
                            "entry for this polytope-fsegment interaction");
                      }
                      ds_segment_prox.erase(finter);
                    },
                    [&]<match::handle<diskfdisk_r> DiskFdiskR>(DiskFdiskR rel) {
                      auto& translat = rel.translated_disk_shape().translation();
                      auto coefs = std::array{translat[0], translat[1]};
                      auto finter = ds_fdisk_prox.find(
                          mp::make_pair(storage::prop<"ds1", Topology>(inter).value(), coefs));
                      if (finter == ds_fdisk_prox.end()) {
                        throw std::runtime_error(
                            "space_filter::remove_interactions: no proximity map "
                            "entry for this ds-fdisk interaction");
                      }
                      ds_fdisk_prox.erase(finter);
                    },
                    [&]<match::handle<polytopefdisk_r<Topology>> PolytopeFdiskR>(
                        PolytopeFdiskR rel) {
                      auto& translat = rel.translated_disk_shape().translation();
                      auto coefs = std::array{translat[0], translat[1]};
                      auto finter = ds_fdisk_prox.find(
                          mp::make_pair(storage::prop<"ds1", Topology>(inter).value(), coefs));
                      if (finter == ds_fdisk_prox.end()) {
                        throw std::runtime_error(
                            "space_filter::remove_interactions: no proximity map "
                            "entry for this polytope-fdisk interaction");
                      }
                      ds_fdisk_prox.erase(finter);
                    },
                    []<bool flag = false>(auto) {
                      // hardening: this used to be a silent no-op in release
                      // (assert compiled out), leaking the interaction
                      throw std::runtime_error(
                          "space_filter::remove_interactions: unhandled relation "
                          "for a self interaction (ds1 == ds2); only "
                          "diskfsegment_r, diskfdisk_r, polytopefsegment_r and "
                          "polytopefdisk_r "
                          "are supported");
                    }));
            //            print("  REMOVE ds segment interaction between {}
            //            {},{},{})\n",
            //                  storage::prop<"ds1", Topology>(inter).get(),
            //                  segmentcoefs[0], segmentcoefs[1],
            //                  segmentcoefs[2], segmentcoefs[3]);

            ;
          }
        }
      }

      auto fact = std::ranges::find(activations, false);

      //      print("  START REMOVE interactions\n");

      while (fact != activations.end()) {
        auto fact_index = fact - activations.begin();
        //        print("  activation of {} is false\n", fact_index);
        auto inter =
            storage::make_handle(data, storage::index<Interaction, indice>(fact_index));

        // with move_back : order is modified

        //        print("  remove interaction {}\n", inter.get());
        storage::remove(data, inter);

        // XXX remove from ds_ds_prox or ds_segment_prox

        // activations has been modified, search first false element
        // starting at current position
        auto remaining_activations = std::ranges::subrange(fact, activations.end());

        auto ifact = std::ranges::find(remaining_activations, false);

        fact = fact + (ifact - remaining_activations.begin());
        // fact = std::ranges::find(activations, false);

        //        print("  find new false at : {}\n", fact -
        //        activations.begin());
      }

      //      print("AFTER REMOVAL: size of indexset0: {}\n",
      //      activations.size()); print("END of interactions removal\n");
      //      print("size of ds ds map: {}\n", ds_ds_prox.size());
      //      print("size of ds segment map: {}\n", ds_segment_prox.size());
    }

    /* remove the ds-mesh interactions the filter did not re-activate this
       step (bodies that left the edge's neighborhood); same move_back
       removal as remove_interactions, the proximity entry goes too */
    void remove_dinteractions(auto& ds_dds_prox) {
      using env = decltype(self()->env());
      using indice = typename env::indice;
      auto& data = self()->data();

      auto& activations = storage::prop_values<dfinteraction, "activation", Topology>(data, 0);
      auto& ds1s = storage::prop_values<dfinteraction, "ds1", Topology>(data, 0);
      auto& ds2s = storage::prop_values<dfinteraction, "ds2", Topology>(data, 0);

      for (auto [activation, ds1, ds2, inter] :
           view::zip(activations, ds1s, ds2s, storage::handles<dfinteraction>(data, 0))) {
        if (activation) continue;
        auto contact_index = variant::visit(data, inter.relation(),
                                            [](auto& rrel) { return rrel.contact_index(); });
        auto finter = ds_dds_prox.find(mp::make_tuple(
            ds1.value(), ds2.value(), edge_point_slot(contact_index, point_of(inter))));
        if (finter == ds_dds_prox.end() ||
            std::get<1>(*finter).value() != inter.index().value()) {
          throw std::runtime_error(
              "space_filter::remove_dinteractions: no proximity map entry for this "
              "ds-mesh interaction (stale filter state or duplicate key)");
        }
        ds_dds_prox.erase(finter);
      }

      // with move_back the last interaction takes the removed index: search
      // again from the same position
      auto fact = std::ranges::find(activations, false);
      while (fact != activations.end()) {
        const auto fact_index = fact - activations.begin();
        storage::remove(data, storage::make_handle(
                                  data, storage::index<dfinteraction, indice>(fact_index)));
        auto remaining = std::ranges::subrange(fact, activations.end());
        fact = fact + (std::ranges::find(remaining, false) - remaining.begin());
      }
    }

    void update_index_set0(auto step) {
      using env = decltype(self()->env());
      using indice = typename env::indice;
      using scalar = typename env::scalar;

      using map_ds_ds_prox_t =
          typename env::template map<mp::tuple<indice, indice, indice>,
                                     storage::index<finteraction, indice>>;

      using map_ds_dds_prox_t =
          typename env::template map<mp::tuple<indice, indice, indice>,
                                     storage::index<dfinteraction, indice>>;

      using map_ds_segment_prox_t =
          typename env::template map<mp::pair<indice, std::array<scalar, 5>>,
                                     storage::index<finteraction, indice>>;

      using map_ds_fdisk_prox_t =
          typename env::template map<mp::pair<indice, std::array<scalar, 2>>,
                                     storage::index<finteraction, indice>>;

      auto& data = self()->data();
      auto& diskfsegments = self()->diskfsegments();

      auto ngbh = self()->neighborhood();
      using ngbh_t = typename std::decay_t<decltype(ngbh)>::type;
      using points_t = typename ngbh_t::points_t;

      auto ds_ds_prox = map_ds_ds_prox_t();
      auto ds_dds_prox = map_ds_dds_prox_t();

      auto ds_segment_prox = map_ds_segment_prox_t();

      auto ds_fdisk_prox = map_ds_fdisk_prox_t();

      const auto n_inter = storage::handles<finteraction>(data, step).size();
      ds_ds_prox.reserve(n_inter * 2);
      ds_segment_prox.reserve(n_inter * 2);
      ds_fdisk_prox.reserve(n_inter * 2);

      auto& ds1s = storage::prop_values<finteraction, "ds1", Topology>(data, step);
      auto& ds2s = storage::prop_values<finteraction, "ds2", Topology>(data, step);

      auto interactions = storage::handles<finteraction>(data, step);

      build_proximity_maps(ds_ds_prox, ds_segment_prox, ds_fdisk_prox, ds1s, ds2s,
                           interactions);

      // no column is attached to empty_item (see storage::item_storage)
      if constexpr (!std::derived_from<dfinteraction, empty_item>) {
        auto& dds1s = storage::prop_values<dfinteraction, "ds1", Topology>(data, step);
        auto& dds2s = storage::prop_values<dfinteraction, "ds2", Topology>(data, step);
        auto dinteractions = storage::handles<dfinteraction>(data, step);
        build_dproximity_maps(ds_dds_prox, dds1s, dds2s, dinteractions);
      }

      auto& activations = storage::prop_values<finteraction, "activation", Topology>(data, 0);

      int activations_size = activations.size();
      if (activations_size > 0) {
        for (auto [activation] : view::zip(activations)) {
          activation = false;
        }
      }

      if constexpr (!std::derived_from<dfinteraction, empty_item>) {
        auto& dactivations =
            storage::prop_values<dfinteraction, "activation", Topology>(data, 0);
        if (dactivations.size() > 0) {
          for (auto [activation] : view::zip(dactivations)) {
            activation = false;
          }
        }
      }

      constexpr auto npointsets = mp::size(points_t{});
      mp::for_each(mp::range<npointsets - mp::size_c<1_c>>, [&](auto ip1) {
        auto p1 = points_t{}[mp::size_c<ip1>];
        using p1_t = decltype(p1);
        auto psid1 = ngbh.point_set_id()[ip1];

        mp::for_each(mp::range_c<std::size_t, ip1, npointsets>, [&](auto ip2) {
          auto p2 = points_t{}[mp::size_c<ip2>];
          using p2_t = decltype(p2);
          auto psid2 = ngbh.point_set_id()[ip2];

          auto& ps1 = ngbh.instance()->point_set(psid1);
          //              auto& ps2 =
          //              ngbh.instance()->point_set(psid2);

          if (ngbh.is_active(ip1, ip2)) {
            for (size_t i = 0; i < ps1.n_points(); ++i) {
              auto pid1 = i;
              auto index_point1 = storage::index<p1_t, size_t>(pid1);
              auto handle_point1 = storage::make_handle(data, index_point1);
              auto body1 = handle_point1.item();

              for (size_t j = 0; j < ps1.n_neighbors(psid2, i); ++j) {
                const unsigned int pid2 = ps1.neighbor(psid2, i, j);

                // print("pid2 : {}\n", pid2);
                auto index_point2 = storage::index<p2_t, size_t>(pid2);
                auto handle_point2 = storage::make_handle(data, index_point2);
                auto body2 = handle_point2.item();
                using system1_t = typename p1_t::item_t;
                using system2_t = typename p2_t::item_t;

                // print("point1 {},{},{}\n",
                //       mp::type_name<system1_t>().c_str(),
                //       handle_point1.coord()(0),
                //       handle_point1.coord()(1));
                // print("point2 {},{},{}\n",
                //       mp::type_name<system2_t>().c_str(),
                //       handle_point2.coord()(0),
                //       handle_point2.coord()(1));

                if constexpr (std::derived_from<system1_t, model::lagrangian_ds>) {
                  // proximity with another disk, only disks are
                  // dynamics check if interaction already exists

                  if constexpr (std::derived_from<system2_t, model::lagrangian_ds>) {
                    dsds_activation(body1, body2, ds_ds_prox, handle_point1, handle_point2);
                  } else {
                    if constexpr (std::derived_from<system2_t, model::elastic_lagrangian_ds>) {
                      dsdds_activation(body1, body2, ds_dds_prox, handle_point1,
                                       handle_point2);
                    } else if constexpr (std::derived_from<system1_t, model::lagrangian_ds>) {
                      if constexpr (std::derived_from<system2_t, collision::shape::segment>) {
                        // body2 is a static segment
                        using body1_t = std::decay_t<decltype(body1)>;
                        if constexpr (collision::always_shape<collision::shape::polytope,
                                                              body1_t>()) {
                          polytopefsegment_activation(body1, body2, ds_segment_prox);
                        } else if constexpr (collision::may_have_shape<
                                                 collision::shape::polytope, body1_t>()) {
                          if (collision::has_shape<collision::shape::polytope>(body1)) {
                            polytopefsegment_activation(body1, body2, ds_segment_prox);
                          } else {
                            diskfsegment_activation(body1, body2, ds_segment_prox);
                          }
                        } else {
                          diskfsegment_activation(body1, body2, ds_segment_prox);
                        }
                      } else if constexpr (std::derived_from<system2_t,
                                                             collision::translated<
                                                                 collision::shape::disk>>) {
                        // body2 is a fixed disk
                        using body1_t = std::decay_t<decltype(body1)>;
                        if constexpr (collision::always_shape<collision::shape::polytope,
                                                              body1_t>()) {
                          polytopefdisk_activation(body1, body2, ds_fdisk_prox);
                        } else if constexpr (collision::may_have_shape<
                                                 collision::shape::polytope, body1_t>()) {
                          if (collision::has_shape<collision::shape::polytope>(body1)) {
                            polytopefdisk_activation(body1, body2, ds_fdisk_prox);
                          } else {
                            diskfdisk_activation(body1, body2, ds_fdisk_prox);
                          }
                        } else {
                          diskfdisk_activation(body1, body2, ds_fdisk_prox);
                        }
                      }
                    }
                  }
                }
              }
            }
          };
        });
      });
      remove_interactions<finteraction>(ds_ds_prox, ds_segment_prox, ds_fdisk_prox);
      if constexpr (!std::derived_from<dfinteraction, empty_item>) {
        remove_dinteractions(ds_dds_prox);
      }
    }

    auto methods() {
      // types only: methods() is called on a bare interface object, without a
      // world (see storage::methods), so nothing here may be evaluated
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      // using scalar = typename env_t::scalar;
      using diskfsegment_r_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<collision::diskfsegment_r, indice>{}))>;
      using diskfdisk_r_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<collision::diskfdisk_r, indice>{}))>;
      using polytopefsegment_r_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(),
          storage::index<collision::polytopefsegment_r<Topology>, indice>{}))>;
      using polytopemesh_r_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(),
          storage::index<collision::polytopemesh_r<Topology>, indice>{}))>;
      using segment_handle_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<collision::shape::segment, indice>{}))>;
      // a handle, not an index: python passes handles, and pybind would
      // reinterpret a handle (data reference first) as an index
      using ref_nslaw_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<nslaw, indice>{}))>;

      return collect(
          method("make_points", &interface<Handle>::make_points),
          method("update_index_set0", &interface<Handle>::update_index_set0<indice>),
          method("insert_diskfsegment_r",
                 &interface<Handle>::insert_diskfsegment_r<diskfsegment_r_t>),
          method("insert_polytopefsegment_r",
                 &interface<Handle>::insert_polytopefsegment_r<polytopefsegment_r_t>),
          method("insert_polytopemesh_r",
                 &interface<Handle>::insert_polytopemesh_r<polytopemesh_r_t>),
          method("insert_diskfdisk_r", &interface<Handle>::insert_diskfdisk_r<diskfdisk_r_t>),
          method("insert_nonsmooth_law",
                 &interface<Handle>::insert_nonsmooth_law<ref_nslaw_t, indice, indice>),
          method("remove_static_segment",
                 &interface<Handle>::template remove_static_item<indice, segment_handle_t>));
    }
  };
};
}  // namespace siconos::collision
