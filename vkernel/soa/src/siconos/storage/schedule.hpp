#pragma once

#include <chrono>
#include <cstddef>
#include <format>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

#include "siconos/storage/kernels.hpp"
#include "siconos/storage/pattern/base.hpp"

namespace siconos::storage::schedule {

/* Schedules: the stages of a computation as a type that reads like the
   algorithm.

   A stage is a struct with a name and a static run(ctx), ctx being the
   context of the computation:

     struct free_state : stage<"free state"> {
       static void run(auto& ctx) { ... }
     };

   A composite has a name, its stages, and a static run(ctx, run_stages)
   that runs them (run_stages(ctx) runs them in order on a context):
   sequence<"name", Stages...> runs them, when<Condition, Stages...> runs
   them if Condition::holds(ctx); an application defines its own (see
   simul::per_element).

   run<S>(ctx, observer) runs a schedule, the observer seeing each stage:
   observer(name, action) must call action() (timers, traces...).
   describe<S>() writes a schedule as text, one stage per line. */

template <pattern::string_literal Name>
struct stage {
  static constexpr std::string_view name{Name.value};
};

// a stage that runs host code besides its kernels (assembled algebra,
// solvers, bookkeeping): the dataflow observer places it on the host
struct host_code {};

// stages in order
template <pattern::string_literal Name, typename... Stages>
struct sequence {
  static constexpr std::string_view name{Name.value};
  using stages = std::tuple<Stages...>;
  static void run(auto& ctx, auto&& run_stages) { run_stages(ctx); }
};

// a condition on the context: a name and a static holds(ctx)
template <pattern::string_literal Name>
struct condition {
  static constexpr std::string_view name{Name.value};
};

// stages run when the condition holds
template <typename Condition, typename... Stages>
struct when {
  static constexpr std::string_view name = Condition::name;
  static constexpr std::string_view prefix = "if ";
  using stages = std::tuple<Stages...>;
  static void run(auto& ctx, auto&& run_stages) {
    if (Condition::holds(ctx)) {
      run_stages(ctx);
    }
  }
};

// an observer that only runs the stages
struct run_only {
  void operator()(std::string_view, auto&& action) const { action(); }
};

// an observer that times the stages: per stage, the number of runs and the
// total time, in the order of their first run. The stages run for each
// element (per_element) are counted together.
struct timings {
  struct entry {
    std::string name;
    std::size_t runs = 0;
    double seconds = 0;
  };
  std::vector<entry> entries;

  void operator()(std::string_view name, auto&& action) {
    const auto start = std::chrono::steady_clock::now();
    action();
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    auto& e = stage(name);
    ++e.runs;
    e.seconds += elapsed.count();
  }

  entry& stage(std::string_view name) {
    for (auto& e : entries) {
      if (e.name == name) return e;
    }
    return entries.emplace_back(entry{std::string(name)});
  }

  double total() const {
    double t = 0;
    for (const auto& e : entries) t += e.seconds;
    return t;
  }

  void clear() { entries.clear(); }

  /// one line per stage: total time, share of the total, runs, name
  std::string report() const {
    const double t = total();
    std::string text =
        std::format("{:>12} {:>7} {:>8}  {}\n", "seconds", "%", "runs", "stage");
    for (const auto& e : entries) {
      text += std::format("{:>12.6f} {:>7.2f} {:>8}  {}\n", e.seconds,
                          t > 0 ? 100. * e.seconds / t : 0., e.runs, e.name);
    }
    text += std::format("{:>12.6f} {:>7.2f} {:>8}  {}\n", t, 100., "", "total");
    return text;
  }
};

// an observer that records what the stages do: the kernels each stage ran,
// in this configuration, with their iteration spaces, the fields they read
// and write and their markers. The report gives a placement for each stage
// (device: its kernels may run on a device; host: a kernel needs host
// objects or an order; no kernel: a host function, the assembled algebra,
// or nothing to do in this configuration) and the candidates for a fusion:
// consecutive kernels on the same items.
struct dataflow {
  struct entry {
    std::string name;
    std::vector<kernels::kernel_info> kernels;
  };
  std::vector<entry> entries;

  void operator()(std::string_view name, auto&& action) {
    std::vector<kernels::kernel_info> ran;
    auto* previous = kernels::recording;
    kernels::recording = &ran;
    action();
    kernels::recording = previous;

    entry* e = nullptr;
    for (auto& x : entries) {
      if (x.name == name) {
        e = &x;
        break;
      }
    }
    if (!e) {
      e = &entries.emplace_back(entry{std::string(name), {}});
    }
    for (auto& k : ran) {
      bool known = false;
      for (const auto& x : e->kernels) {
        known = known || (x.name == k.name && x.item == k.item);
      }
      if (!known) e->kernels.push_back(std::move(k));
    }
  }

  static std::string placement(const entry& e) {
    if (e.kernels.empty()) return "no kernel";
    for (const auto& k : e.kernels) {
      if (k.host_only || k.ordered) return "host";
    }
    return "device";
  }

  static bool fusable(const kernels::kernel_info& a, const kernels::kernel_info& b) {
    return a.space == b.space && a.item == b.item && !a.reads_world && !b.reads_world &&
           !a.ordered && !b.ordered;
  }

  std::string report() const {
    auto join = [](const std::vector<std::string>& v) {
      std::string t;
      for (const auto& x : v) t += (t.empty() ? "" : ", ") + x;
      return t;
    };
    std::string text;
    for (const auto& e : entries) {
      text += std::format("[{}] {}\n", placement(e), e.name);
      for (const auto& k : e.kernels) {
        text += std::format("    {}: {} ({})\n", k.name, k.space, k.item);
        text += std::format("      reads {}\n", k.reads.empty() ? "-" : join(k.reads));
        auto writes = k.writes;
        if (k.assembles) writes.push_back("the assembled problem (by rank)");
        text += std::format("      writes {}\n", writes.empty() ? "-" : join(writes));
        if (k.host_only || k.ordered || k.reads_world) {
          text += std::format("      {}{}{}\n", k.host_only ? "host only " : "",
                              k.ordered ? "ordered " : "", k.reads_world ? "reads world" : "");
        }
      }
    }
    text += "fusion candidates:\n";
    const kernels::kernel_info* previous = nullptr;
    for (const auto& e : entries) {
      // a stage without kernel (a host function) is a barrier
      if (e.kernels.empty()) previous = nullptr;
      for (const auto& k : e.kernels) {
        if (previous && fusable(*previous, k)) {
          text += std::format("  {} | {}: {} ({})\n", previous->name, k.name, k.space, k.item);
        }
        previous = &k;
      }
    }
    return text;
  }
};

// the timings and the dataflow, each one optional (null): one observer type
// whatever is enabled, the schedule is instantiated once
struct optional_observers {
  timings* times = nullptr;
  dataflow* flow = nullptr;
  void operator()(std::string_view name, auto&& action) {
    if (times && flow) {
      (*times)(name, [&]() { (*flow)(name, action); });
    } else if (times) {
      (*times)(name, action);
    } else if (flow) {
      (*flow)(name, action);
    } else {
      action();
    }
  }
};

// two observers: a sees the stage, then b
template <typename A, typename B>
struct both {
  A& a;
  B& b;
  void operator()(std::string_view name, auto&& action) {
    a(name, [&]() { b(name, action); });
  }
};

template <typename S>
static constexpr bool is_composite = requires { typename S::stages; };

/// run the schedule S on the context ctx
template <typename S>
void run(auto& ctx, auto&& observer) {
  if constexpr (is_composite<S>) {
    S::run(ctx, [&observer]<typename Ctx>(Ctx& c) {
      [&]<typename... Stages>(std::tuple<Stages...>*) {
        (run<Stages>(c, observer), ...);
      }(static_cast<typename S::stages*>(nullptr));
    });
  } else {
    observer(S::name, [&ctx]() {
      if constexpr (std::is_base_of_v<host_code, S>) {
        if (kernels::recording) {
          kernels::recording->push_back(kernels::kernel_info{
              .name = "host code", .space = "-", .item = "-", .host_only = true});
        }
      }
      S::run(ctx);
    });
  }
}

template <typename S>
void run(auto& ctx) {
  run<S>(ctx, run_only{});
}

/// the schedule S as text, one stage per line, the stages of a composite
/// indented under its name
template <typename S>
std::string describe(std::size_t depth = 0) {
  std::string text(2 * depth, ' ');
  if constexpr (is_composite<S>) {
    if constexpr (requires { S::prefix; }) {
      text += S::prefix;
    }
    text += S::name;
    text += ":\n";
    [&]<typename... Stages>(std::tuple<Stages...>*) {
      ((text += describe<Stages>(depth + 1)), ...);
    }(static_cast<typename S::stages*>(nullptr));
  } else {
    text += S::name;
    text += "\n";
  }
  return text;
}

}  // namespace siconos::storage::schedule
