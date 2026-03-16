#include <numeric>
#include <print>
#include <random>

#include "ProblemMatrix.h"
#include "seekers/BruteForce.h"
#include "seekers/SimHash.h"
#include "seekers/Tables.h"

static const std::vector<std::string> problems = {"markshare_4_0",
                                                  "markshare2",
                                                  "gen-ip054",
                                                  "enlight_hard",
                                                  "timtab1",
                                                  "pk1",
                                                  "gen-ip002",
                                                  "sp150x300d",
                                                  "neos-3046615-murg",
                                                  "neos859080",
                                                  "neos-3754480-nidda",
                                                  "mas76",
                                                  "neos-2657525-crna",
                                                  "mas74",
                                                  "glass4",
                                                  "exp-1-500-5-5",
                                                  "neos5",
                                                  "fhnw-binpack4-4",
                                                  "supportcase26",
                                                  "tr12-30",
                                                  "graphdraw-domain",
                                                  "50v-10",
                                                  "mad",
                                                  "ic97_potential",
                                                  "neos-911970",
                                                  "assign1-5-8",
                                                  "neos-3381206-awhea",
                                                  "binkar10_1",
                                                  "p200x1188c",
                                                  "gmu-35-40",
                                                  "neos17",
                                                  "beasleyC3",
                                                  "neos-3627168-kasai",
                                                  "pg",
                                                  "csched008",
                                                  "mc11",
                                                  "neos-4338804-snowy",
                                                  "csched007",
                                                  "lotsize",
                                                  "pg5_34",
                                                  "neos-4954672-berkel",
                                                  "mcsched",
                                                  "gmu-35-50",
                                                  "mik-250-20-75-4",
                                                  "n5-3",
                                                  "uct-subprob",
                                                  "ran14x18-disj-8",
                                                  "b1c1s1",
                                                  "rococoC10-001000",
                                                  "cost266-UUE",
                                                  "cvs16r128-89",
                                                  "reblock115",
                                                  "milo-v12-6-r2-40-1",
                                                  "rococoB10-011000",
                                                  "fhnw-binpack4-48",
                                                  "neos-3024952-loue",
                                                  "nu25-pr12",
                                                  "qap10",
                                                  "app1-1",
                                                  "graph20-20-1rand",
                                                  "fastxgemm-n2r6s0t2",
                                                  "neos-5107597-kakapo",
                                                  "neos-4738912-atrato",
                                                  "bnatt400",
                                                  "rmatr100-p10",
                                                  "icir97_tension",
                                                  "neos-4387871-tavua",
                                                  "sct2",
                                                  "bppc4-08",
                                                  "neos-3083819-nubu",
                                                  "neos-1582420",
                                                  "bnatt500",
                                                  "rocI-4-11",
                                                  "peg-solitaire-a3",
                                                  "supportcase18",
                                                  "roll3000",
                                                  "h80x6320d",
                                                  "traininstance6",
                                                  "seymour",
                                                  "seymour1",
                                                  "swath1",
                                                  "swath3",
                                                  "neos-1456979",
                                                  "neos-1445765",
                                                  "traininstance2",
                                                  "drayage-100-23",
                                                  "drayage-25-23",
                                                  "neos-2978193-inde",
                                                  "neos-3004026-krka",
                                                  "nexp-150-20-8-5",
                                                  "eil33-2",
                                                  "CMS750_4",
                                                  "piperout-08",
                                                  "air05",
                                                  "piperout-27",
                                                  "dws008-01",
                                                  "comp21-2idx",
                                                  "cod105",
                                                  "neos-1171737",
                                                  "neos-3656078-kumeu",
                                                  "glass-sc",
                                                  "decomp2",
                                                  "cbs-cta",
                                                  "istanbul-no-cutoff",
                                                  "dano3_3",
                                                  "dano3_5",
                                                  "net12",
                                                  "ns1208400",
                                                  "comp07-2idx",
                                                  "wachplan",
                                                  "neos-3216931-puriri",
                                                  "trento1",
                                                  "radiationm18-12-05",
                                                  "irp",
                                                  "ns1830653",
                                                  "momentum1",
                                                  "neos-950242",
                                                  "supportcase40",
                                                  "fiball"};

int main() {
  std::string problem_name = "neos-5114902-kasavu";
  // for (const auto& problem_name : problems) {
  auto matrix = get_problem_matrix(problem_name);

  std::println("{}: {} x {}", problem_name, matrix.shape().first,
               matrix.shape().second);

  // auto bf_result =
  // seekers::BruteForce(seekers::BruteForceHamming{2}).seek(matrix);
  auto result = seekers::Tables().seek(matrix);
  std::println("found: {} pairs", result.size());
  // }

  // for (auto [i, j] : result) {
  //   std::println("------------- ({}, {}) -------------", i, j);
  //
  //   auto i_row = matrix.get_row(i);
  //   auto j_row = matrix.get_row(j);
  //
  //   for (auto [col, value] : i_row) {
  //     std::print("({}, {}) ", col, value);
  //   }
  //   std::print("\n");
  //
  //   for (auto [col, value] : j_row) {
  //     std::print("({}, {}) ", col, value);
  //   }
  //   std::print("\n\n");
  // }

  // auto sh_result = seekers::SimHash().seek(matrix, bf_result);
  // std::println("found: {} pairs", sh_result.size());

  // for (const auto& problem_name : problems) {
  //   auto matrix = get_problem_matrix(problem_name);
  //
  //   auto [n, d] = matrix.shape();
  //
  //   std::vector<size_t> permutation(d);
  //   std::iota(permutation.begin(), permutation.end(), 0);
  //   std::ranges::shuffle(permutation, std::default_random_engine{});
  //
  //   std::vector<size_t> counts(n);
  //   for (size_t i = 0; 2 * i < d; ++i) {
  //     for (auto [row, _] : matrix.get_column(permutation[i])) {
  //       ++counts[row];
  //     }
  //   }
  //
  //   size_t zero_count = 0;
  //   for (size_t cnt : counts) {
  //     if (cnt == 0) {
  //       ++zero_count;
  //     }
  //   }
  //
  //   std::println("{} - {}/{} = {}", problem_name, zero_count, n,
  //   static_cast<double>(zero_count) / static_cast<double>(n));
  // }
}
