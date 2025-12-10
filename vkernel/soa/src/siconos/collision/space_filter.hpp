#pragma once

#include <concepts>
#include <range/v3/view/enumerate.hpp>

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/diskdisk_r.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/diskfsegment_r.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/chained_segment.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/variant.hpp"

namespace std {
template <typename T>
struct hash<std::array<T, 2>> {
  size_t operator()(const std::array<T, 2>& arr) const noexcept
  {
    size_t h1 = std::hash<T>{}(arr[0]);
    size_t h2 = std::hash<T>{}(arr[1]);
    return h1 ^ (h2 << 1);
  }
};

template <typename T>
struct hash<std::array<T, 4>> {
  size_t operator()(const std::array<T, 4>& arr) const noexcept
  {
    size_t h1 = std::hash<T>{}(arr[0]);
    size_t h2 = std::hash<T>{}(arr[1]);
    size_t h3 = std::hash<T>{}(arr[2]);
    size_t h4 = std::hash<T>{}(arr[3]);
    size_t combined1 = h1 ^ (h2 << 1);
    size_t combined2 = h3 ^ (h4 << 1);
    return combined1 ^ (combined2 << 1);
  }
};

template <typename First, typename Second>
struct hash<siconos::storage::ground::pair<First, Second>> {
  size_t operator()(
      const siconos::storage::ground::pair<First, Second>& p) const noexcept
  {
    size_t h1 = std::hash<First>{}(siconos::storage::ground::first(p));
    size_t h2 = std::hash<Second>{}(siconos::storage::ground::second(p));
    return h1 ^ (h2 << 1);
  }
};
}  // namespace std

namespace siconos::collision {

template <typename Topology, typename Neighborhood>
struct space_filter : item<> {
  using items = gather<Topology, Neighborhood>;
  using topology = Topology;
  using dynamical_system = typename topology::fixed_dof_system;
  using interaction = typename topology::finteraction;
  using nslaw = typename interaction::nslaw;

  using neighborhood = Neighborhood;

  using attributes = gather<
      attribute<"topology", some::item_ref<topology>>,
      attribute<"neighborhood", some::item_ref<neighborhood>>,
      attribute<"nslaw", some::item_ref<nslaw>>,
      attribute<"diskdisk_r", some::item_ref<diskdisk_r>>,
      attribute<"diskmeshes",
                some::map<some::array<some::indice, some::indice_value<2>>,
                          some::item_ref<diskmesh_r>>>,
      attribute<"diskfsegments",
                some::map<some::array<some::scalar, some::indice_value<4>>,
                          some::item_ref<diskfsegment_r>>>,
      attribute<"diskfdisks",
                some::map<some::array<some::scalar, some::indice_value<2>>,
                          some::item_ref<diskfdisk_r>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nslaw() { return storage::attr<"nslaw">(*self()); }
    decltype(auto) topology()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"topology">(*self()));
    }

    decltype(auto) neighborhood()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"neighborhood">(*self()));
    }

    decltype(auto) diskdisk_r()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"diskdisk_r">(*self()));
    }

    decltype(auto) diskmeshes()
    {
      return storage::attr<"diskmeshes">(*self());
    }

    decltype(auto) diskfsegments()
    {
      return storage::attr<"diskfsegments">(*self());
    }

    decltype(auto) diskfdisks()
    {
      return storage::attr<"diskfdisks">(*self());
    }

    // void create_relations()
    // {
    //   using env_t = decltype(self()->env());
    //   using indice = typename env_t::indice;

    //   using points_t = typename decltype(self()->neighborhood())::points_t;

    //   auto& data = self()->data();

    //   ground::for_each(points_t{}, []<typename Point>(Point) {
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

    void remove_static_item(auto step, auto item_handle)
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using item_t = typename std::decay_t<decltype(item_handle)>::type;
      using point_t = collision::point<item_t>;
      using points_t = typename neighborhood::points_t;

      auto& data = self()->data();

      auto& points_flags = storage::attr_values<point_t, "flag">(data, step);
      auto& points_items = storage::attr_values<point_t, "item">(data, step);
      auto& points_coords =
          storage::attr_values<point_t, "coord">(data, step);

      auto ps_indx = ground::index_of<point_t>(ground::std_tuple(points_t{}));

      // ground::index_if(
      //   points_t{},
      //   ground::equal.to(collision::point<item_t>));

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
        auto point = storage::make_handle(
            data, storage::index<point_t, indice>(ff_index));
        storage::remove(data, point);

        auto remaining_points_flags =
            std::ranges::subrange(ff, points_flags.end());
        auto rff = std::ranges::find(remaining_points_flags, true);

        ff = ff + (rff - remaining_points_flags.begin());
      };

      neighborhood().instance()->resize_point_set(
          ps_indx, points_coords.front().data(), points_coords.size());
    }

    void insert_diskfsegment_r(auto hdl)
    {
      storage::attr<"diskfsegments">(
          *self())[{hdl.segment().x1(), hdl.segment().x2(),
                    hdl.segment().y1(), hdl.segment().y2()}] = hdl;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void insert_diskfdisk_r(auto tds)
    {
      auto htds = tds.translated_disk_shape();
      storage::attr<"diskfdisks">(
          *self())[{htds.translation()[0], htds.translation()[1]}] = tds;

      // for (auto hds : storage::handles<dynamical_system>(self()->data())) {
      //   auto inter = self()->topology().link(hds);
      //   inter.relation() = dl;
      //   inter.nslaw() = nslaw();
      // }
    }

    void make_points()
    {
      auto& data = self()->data();
      using env = decltype(self()->env());
      using indice = typename env::indice;

      using points_t = typename neighborhood::points_t;

      ground::for_each(points_t{}, [&]<typename Point>(Point) {
        using item_t = typename Point::item_t;

        if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
          auto all_ds = storage::handles<item_t>(data);
          for (auto ds : all_ds) {
            //            print("add disk point : {}\n", ds.get());
            auto shape = storage::prop<"shape">(ds);
            using shape_t = std::decay_t<decltype(shape)>;
            if constexpr (std::derived_from<
                              shape_t, collision::shape::chained_segment>) {
              auto nbsegments = std::size(shape.points()) / 2;
              // multiple points are associated to the system
              for (std::size_t index = 0; index < nbsegments; ++index) {
                for (auto [i, point_coord] :
                     shape.points_coords(index) | view::enumerate) {
                  auto new_point = storage::add<Point>(data);
                  new_point.item() = ds;
                  new_point.coord() = point_coord;
                  new_point.seg_index() = index;
                  new_point.point_index() = i;
                }
              }
            }
            else {
              // the system has just one point associated
              auto new_point = storage::add<Point>(data);
              new_point.item() = ds;
              new_point.update(0);
            }
          }
        }
        // fixed segment
        else if constexpr (std::derived_from<item_t,
                                             collision::shape::segment>) {
          auto all_segments = storage::handles<item_t>(data);
          for (auto segment : all_segments) {
            for (auto [i, point_coord] :
                 segment.points_coords() | view::enumerate) {
              auto new_point = storage::add<Point>(data);
              new_point.item() = segment;
              new_point.coord() = point_coord;
              new_point.point_index() = i;
            }
          }
        }
        // fixed disk
        else if constexpr (std::derived_from<item_t,
                                             collision::translated<
                                                 collision::shape::disk>>) {
          auto all_fdisks = storage::handles<item_t>(data);
          for (auto fdisk : all_fdisks) {
            for (auto [i, point_coord] :
                 fdisk.points_coords() | view::enumerate) {
              auto new_point = storage::add<Point>(data);
              new_point.item() = fdisk;
              new_point.coord() = point_coord;
              new_point.point_index() = (indice) i;
            }
          }
        };
      });
    }

    decltype(auto) make_ipair(auto ids1, auto ids2)
    {
      auto i1 = ids1.value();
      auto i2 = ids2.value();
      if (i1 < i2) {
        return ground::make_pair(i1, i2);
      }
      else {
        return ground::make_pair(i2, i1);
      }
    }

    void build_proximity_maps(auto& ds_ds_prox, auto& ds_segment_prox,
                              auto& ds_fdisk_prox, const auto& ds1s,
                              const auto& ds2s, const auto& interactions)
    {
      auto& data = self()->data();

      // build (ds ds -> inter) & (ds segment -> inter) maps
      for (auto [ds1, ds2, inter] : view::zip(ds1s, ds2s, interactions)) {
        if (ds1 != ds2) {
          ds_ds_prox[make_ipair(ds1, ds2)] = inter;
        }
        else {
          siconos::variant::visit(
              data, inter.relation(),
              ground::overload(
                  // https://stackoverflow.com/questions/46114214/lambda-implicit-capture-fails-with-variable-declared-from-structured-binding
                  // capture by value ok with handles
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<diskfsegment_r> DiskSegmentR>(
                      DiskSegmentR rel) {
                    auto segment = rel.segment();
                    auto coefs = std::array{segment.x1(), segment.x2(),
                                            segment.y1(), segment.y2()};
                    ds_segment_prox[ground::make_pair(ds1.value(), coefs)] =
                        inter;
                  },
                  [&, ds1 = ds1,
                   inter = inter]<match::handle<diskfdisk_r> DiskFdiskR>(
                      DiskFdiskR rel) {
                    auto fdisk = rel.translated_disk_shape();
                    auto& translat = fdisk.translation();
                    auto coefs = std::array{translat[0], translat[1]};
                    ds_fdisk_prox[ground::make_pair(ds1.value(), coefs)] =
                        inter;
                  },
                  []<bool flag = false>(auto) {
                    assert(flag);
                    // should send an exception here
                    // static_assert(flag,
                    //               "should not
                    //               happen");
                  }));
        }
      }
    }

    void dsds_activation(auto& ds1, auto& ds2, auto& map, auto hp1, auto hp2)
    {
      auto& data = self()->data();
      auto topo = self()->topology();
      auto nslaw = self()->nslaw();
      auto diskdisk_r = self()->diskdisk_r();
      auto relmap = self()->diskmeshes();

      // at most one edge between 2 ds !!
      auto find_inter = map.find(make_ipair(ds1.index(), ds2.index()));
      if (find_inter != map.end()) {
        // keep this edge
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        storage::prop<"activation">(inter) = true;
      }
      else {
        // create the edge
        auto inter = topo.link(ds1, ds2);
        inter.nslaw() = nslaw;  // one nslaw for the moment

        storage::prop<"activation">(inter) = true;

        auto shape1 = storage::prop<"shape">(ds1);
        using shape1_t = std::decay_t<decltype(shape1)>;

        auto shape2 = storage::prop<"shape">(ds2);
        using shape2_t = std::decay_t<decltype(shape2)>;

        if constexpr (std::derived_from<shape1_t, collision::shape::disk>) {
          if constexpr (std::derived_from<
                            shape2_t, collision::shape::chained_segment>) {
            auto mesh = hp2.mesh();
            auto contact_index = hp2.seg_index();
            if (auto search = relmap.find({mesh.get(), contact_index});
                search != relmap.end()) {
              inter.relation = search->second;
            }
            else {
              auto rel = storage::add<collision::diskmesh_r>(data);
              rel.mesh() = mesh;
              rel.contact_index = contact_index;
              inter.relation() = rel;
              relmap[{mesh.get(), contact_index}] = rel;
            }
          }
          else if constexpr (std::derived_from<shape2_t,
                                               collision::shape::disk>) {
            inter.relation() = diskdisk_r;
          }
          else {
            assert(false);
          }
        }

        map[make_ipair(ds1.index(), ds2.index())] = inter;
      }
    }

    void diskfsegment_activation(auto& body1, auto& body2, auto& intermap)
    {
      auto& data = self()->data();
      auto nslaw = self()->nslaw();
      auto topo = self()->topology();
      auto& relmap = self()->diskfsegments();

      auto segment = body2;
      auto x1 = segment.x1();
      auto x2 = segment.x2();
      auto y1 = segment.y1();
      auto y2 = segment.y2();

      auto& ds1 = body1;

      auto find_inter = intermap.find(
          ground::make_pair(ds1.index().value(), std::array{x1, x2, y1, y2}));

      if (find_inter != intermap.end()) {
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        // keep this interaction
        storage::prop<"activation">(inter) = true;

        // print("interaction FOUND for
        // {},{},{}\n", a, b,
        //       c);
      }
      else {
        // print("interaction NOT FOUND for
        // {},{},{}\n", a,
        //       b, c);
        // create the edge
        auto inter = topo.link(body1);
        inter.nslaw() = nslaw;  // one nslaw for the moment

        if (auto search = relmap.find({x1, x2, y1, y2});
            search != relmap.end()) {
          inter.relation() = search->second;
        }
        else {
          auto dl = storage::add<diskfsegment_r>(data);
          auto segment = dl.segment();
          segment.x1() = x1;
          segment.x2() = x2;
          segment.y1() = y1;
          segment.y2() = y2;

          inter.relation() = dl;
          relmap[{x1, x2, y1, y2}] = dl;
        }
        storage::prop<"activation">(inter) = true;
        intermap[ground::make_pair(ds1.index().value(),
                                   std::array{x1, x2, y1, y2})] = inter;
      }
    }

    void diskfdisk_activation(auto& body1, auto& body2, auto& ds_fdisk_prox)
    {
      auto& data = self()->data();
      auto nslaw = self()->nslaw();
      auto topo = self()->topology();
      auto diskfdisks = self()->diskfdisks();

      auto fdisk = body2;
      auto& translat = fdisk.translation();
      auto coefs = std::array{translat[0], translat[1]};

      auto& ds1 = body1;

      auto find_inter =
          ds_fdisk_prox.find(ground::make_pair(ds1.index().value(), coefs));

      if (find_inter != ds_fdisk_prox.end()) {
        auto inter = storage::make_handle(data, std::get<1>(*find_inter));
        storage::prop<"activation">(inter) = true;
      }
      else {
        auto inter = topo.link(body1);
        inter.nslaw() = nslaw;

        if (auto search = diskfdisks.find(coefs);
            search != diskfdisks.end()) {
          inter.relation() = search->second;
        }
        else {
          auto dfd = storage::add<diskfdisk_r>(data);
          auto tds = dfd.translated_disk_shape();
          tds.translation() = translat;

          inter.relation() = dfd;
          diskfdisks[coefs] = dfd;
        }
        storage::prop<"activation">(inter) = true;
        ds_fdisk_prox[ground::make_pair(ds1.index().value(), coefs)] = inter;
      }
    }

    void remove_interactions(auto& ds_ds_prox, auto& ds_segment_prox,
                             auto& ds_fdisk_prox)
    {
      using env = decltype(self()->env());
      using indice = typename env::indice;

      auto& data = self()->data();

      auto& activations =
          storage::prop_values<interaction, "activation">(data, 0);

      auto interactions = storage::handles<interaction>(data, 0);

      //      print("BEFORE REMOVAL: size of indexset0: {}\n",
      //      activations.size());
      for (auto [activation, inter] : view::zip(activations, interactions)) {
        if (!activation) {
          //          print("START REMOVE interaction {}\n", inter.get());

          if (storage::prop<"ds1">(inter) != storage::prop<"ds2">(inter)) {
            auto finter = ds_ds_prox.find(make_ipair(
                storage::prop<"ds1">(inter), storage::prop<"ds2">(inter)));

            assert(inter.index().value() == std::get<1>(*finter).value());

            ds_ds_prox.erase(finter);
            //            print("  REMOVE ds ds interaction between {}
            //            {}\n",
            //                  storage::prop<"ds1">(inter).get(),
            //                  storage::prop<"ds2">(inter).get());
          }
          else {
            siconos::variant::visit(
                data, inter.relation(),
                ground::overload(
                    [&]<match::handle<diskfsegment_r> DiskSegmentR>(
                        DiskSegmentR rel) {
                      auto segment = rel.segment();
                      auto coefs = std::array{segment.x1(), segment.x2(),
                                              segment.y1(), segment.y2()};
                      auto finter = ds_segment_prox.find(ground::make_pair(
                          storage::prop<"ds1">(inter).value(), coefs));
                      ds_segment_prox.erase(finter);
                    },
                    [&]<match::handle<diskfdisk_r> DiskFdiskR>(
                        DiskFdiskR rel) {
                      auto& translat =
                          rel.translated_disk_shape().translation();
                      auto coefs = std::array{translat[0], translat[1]};
                      auto finter = ds_fdisk_prox.find(ground::make_pair(
                          storage::prop<"ds1">(inter).value(), coefs));
                      ds_fdisk_prox.erase(finter);
                    },
                    []<bool flag = false>(auto) {
                      assert(flag);
                      // should send an exception here
                      // static_assert(flag,
                      //               "should not
                      //               happen");
                    }));
            //            print("  REMOVE ds segment interaction between {}
            //            {},{},{})\n",
            //                  storage::prop<"ds1">(inter).get(),
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
        auto inter = storage::make_handle(
            data, storage::index<interaction, indice>(fact_index));

        // with move_back : order is modified

        //        print("  remove interaction {}\n", inter.get());
        storage::remove(data, inter);

        // XXX remove from ds_ds_prox or ds_segment_prox

        // activations has been modified, search first false element
        // starting at current position
        auto remaining_activations =
            std::ranges::subrange(fact, activations.end());

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

    void update_index_set0(auto step)
    {
      using env = decltype(self()->env());
      using indice = typename env::indice;
      using scalar = typename env::scalar;

      using map_ds_ds_prox_t =
          typename env::template map<ground::pair<indice, indice>,
                                     storage::index<interaction, indice>>;

      using map_ds_segment_prox_t = typename env::template map<
          ground::pair<indice, std::array<scalar, 4>>,
          storage::index<interaction, indice>>;

      using map_ds_fdisk_prox_t = typename env::template map<
          ground::pair<indice, std::array<scalar, 2>>,
          storage::index<interaction, indice>>;

      auto& data = self()->data();
      auto diskfsegments = self()->diskfsegments();

      auto ngbh = self()->neighborhood();
      using ngbh_t = typename std::decay_t<decltype(ngbh)>::type;
      using points_t = typename ngbh_t::points_t;

      auto ds_ds_prox = map_ds_ds_prox_t();

      auto ds_segment_prox = map_ds_segment_prox_t();

      auto ds_fdisk_prox = map_ds_fdisk_prox_t();

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto interactions = storage::handles<interaction>(data, step);

      build_proximity_maps(ds_ds_prox, ds_segment_prox, ds_fdisk_prox, ds1s,
                           ds2s, interactions);

      auto& activations =
          storage::prop_values<interaction, "activation">(data, 0);

      int activations_size = activations.size();
      if (activations_size > 0) {
        for (auto [activation] : view::zip(activations)) {
          activation = false;
        }
      }

      constexpr auto npointsets = ground::size(points_t{});
      ground::for_each(
          ground::range<npointsets - ground::size_c<1_c>>, [&](auto ip1) {
            auto p1 = points_t{}[ground::size_c<ip1>];
            using p1_t = decltype(p1);
            auto psid1 = ngbh.point_set_id()[ip1];

            ground::for_each(
                ground::range_c<std::size_t, ip1, npointsets>, [&](auto ip2) {
                  auto p2 = points_t{}[ground::size_c<ip2>];
                  using p2_t = decltype(p2);
                  auto psid2 = ngbh.point_set_id()[ip2];

                  auto& ps1 = ngbh.instance()->point_set(psid1);
                  //              auto& ps2 =
                  //              ngbh.instance()->point_set(psid2);

                  if (ngbh.is_active(ip1, ip2)) {
                    for (size_t i = 0; i < ps1.n_points(); ++i) {
                      auto pid1 = i;
                      auto index_point1 = storage::index<p1_t, size_t>(pid1);
                      auto handle_point1 =
                          storage::make_handle(data, index_point1);
                      auto body1 = handle_point1.item();

                      for (size_t j = 0; j < ps1.n_neighbors(psid2, i); ++j) {
                        const unsigned int pid2 = ps1.neighbor(psid2, i, j);

                        // print("pid2 : {}\n", pid2);
                        auto index_point2 =
                            storage::index<p2_t, size_t>(pid2);
                        auto handle_point2 =
                            storage::make_handle(data, index_point2);
                        auto body2 = handle_point2.item();
                        using system1_t = typename p1_t::item_t;
                        using system2_t = typename p2_t::item_t;

                        // print("point1 {},{},{}\n",
                        //       ground::type_name<system1_t>().c_str(),
                        //       handle_point1.coord()(0),
                        //       handle_point1.coord()(1));
                        // print("point2 {},{},{}\n",
                        //       ground::type_name<system2_t>().c_str(),
                        //       handle_point2.coord()(0),
                        //       handle_point2.coord()(1));

                        if constexpr (std::derived_from<
                                          system1_t, model::lagrangian_ds>) {
                          // proximity with another disk, only disks are
                          // dynamics check if interaction already exists

                          if constexpr (std::derived_from<
                                            system2_t,
                                            model::lagrangian_ds>) {
                            dsds_activation(body1, body2, ds_ds_prox,
                                            handle_point1, handle_point2);
                          }
                          else {
                            if constexpr (std::derived_from<
                                              system1_t,
                                              model::lagrangian_ds>) {
                              if constexpr (std::derived_from<
                                                system2_t,
                                                collision::shape::segment>) {
                                // body2 is a static segment
                                // for all self edges find the one with the
                                // corresponding segment
                                diskfsegment_activation(body1, body2,
                                                        ds_segment_prox);
                              }
                              else if constexpr (std::derived_from<
                                                     system2_t,
                                                     collision::translated<
                                                         collision::shape::
                                                             disk>>) {
                                diskfdisk_activation(body1, body2,
                                                     ds_fdisk_prox);
                              }
                            }
                          }
                        }
                      }
                    }
                  };
                });
          });
      remove_interactions(ds_ds_prox, ds_segment_prox, ds_fdisk_prox);
    }

    auto methods()
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      // using scalar = typename env_t::scalar;
      using diskfsegment_r_t = std::decay_t<decltype(storage::make_handle(
          data, storage::index<collision::diskfsegment_r, indice>{}))>;
      using diskfdisk_r_t = std::decay_t<decltype(storage::make_handle(
          data, storage::index<collision::diskfdisk_r, indice>{}))>;
      using segment_handle_t = std::decay_t<decltype(storage::make_handle(
          data, storage::index<collision::shape::segment, indice>{}))>;

      return collect(
          method("make_points", &interface<Handle>::make_points),
          method("update_index_set0",
                 &interface<Handle>::update_index_set0<indice>),
          method("insert_diskfsegment_r",
                 &interface<Handle>::insert_diskfsegment_r<diskfsegment_r_t>),
          method("insert_diskfdisk_r",
                 &interface<Handle>::insert_diskfdisk_r<diskfdisk_r_t>),
          method("remove_static_segment",
                 &interface<Handle>::template remove_static_item<
                     indice, segment_handle_t>));
    }
  };
};
}  // namespace siconos::collision
