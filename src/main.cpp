#include <numeric>
#include <print>
#include <random>

#include "ProblemMatrix.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "similarity/ZipRows.h"

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
                                                  "fiball",
                                                  "30n20b8",
                                                  "neos-827175",
                                                  "rmatr200-p5",
                                                  "satellites2-60-fs",
                                                  "unitcal_7",
                                                  "leo1",
                                                  "neos-1171448",
                                                  "brazil3",
                                                  "neos-5188808-nattai",
                                                  "mzzv11",
                                                  "chromaticindex512-7",
                                                  "mzzv42z",
                                                  "cmflsp50-24-8-8",
                                                  "neos-1122047",
                                                  "germanrr",
                                                  "neos-5195221-niemur",
                                                  "neos-933966",
                                                  "n2seq36q",
                                                  "neos-4300652-rahue",
                                                  "neos-1354092",
                                                  "mushroom-best",
                                                  "neos-960392",
                                                  "blp-ic98",
                                                  "app1-2",
                                                  "neos-662469",
                                                  "blp-ar98",
                                                  "nursesched-sprint02",
                                                  "supportcase33",
                                                  "thor50dday",
                                                  "tbfp-network",
                                                  "leo2",
                                                  "lectsched-5-obj",
                                                  "atlanta-ip",
                                                  "neos-4763324-toguru",
                                                  "sing326",
                                                  "academictimetablesmall",
                                                  "chromaticindex1024-7",
                                                  "neos-2746589-doon",
                                                  "sing44",
                                                  "satellites2-40",
                                                  "sp97ar",
                                                  "cryptanalysiskb128n5obj14",
                                                  "cryptanalysiskb128n5obj16",
                                                  "neos-787933",
                                                  "rocII-5-11",
                                                  "neos-4722843-widden",
                                                  "neos8",
                                                  "uccase9",
                                                  "ns1952667",
                                                  "sorrell3",
                                                  "n3div36",
                                                  "neos-3988577-wolgan",
                                                  "neos-873061",
                                                  "opm2-z10-s4",
                                                  "neos-860300",
                                                  "rail01",
                                                  "radiationm40-10-02",
                                                  "fast0507",
                                                  "ns1116954",
                                                  "uccase12",
                                                  "sp98ar",
                                                  "hypothyroid-k1",
                                                  "supportcase42",
                                                  "snp-02-004-104",
                                                  "rail507",
                                                  "physiciansched6-2",
                                                  "neos-957323",
                                                  "proteindesign122trx11p8",
                                                  "triptim1",
                                                  "ex9",
                                                  "irish-electricity",
                                                  "neos-4532248-waihi",
                                                  "map10",
                                                  "map16715-04",
                                                  "supportcase10",
                                                  "neos-2987310-joes",
                                                  "supportcase6",
                                                  "netdiversion",
                                                  "nursesched-medium-hint03",
                                                  "proteindesign121hz512p9",
                                                  "nw04",
                                                  "rail02",
                                                  "neos-4413714-turia",
                                                  "neos-5093327-huahum",
                                                  "buildingenergy",
                                                  "neos-3555904-turama",
                                                  "neos-631710",
                                                  "var-smallemery-m6j6",
                                                  "rd-rplusc-21",
                                                  "roi2alpha3n4",
                                                  "eilA101-2",
                                                  "physiciansched3-3",
                                                  "neos-848589",
                                                  "ex10",
                                                  "gfd-schedulen180f7d50m30k18",
                                                  "bab6",
                                                  "s250r10",
                                                  "neos-5049753-cuanza",
                                                  "highschool1-aigio",
                                                  "k1mushroom",
                                                  "splice1k1",
                                                  "savsched1",
                                                  "s100",
                                                  "ns1760995",
                                                  "co-100",
                                                  "bab2",
                                                  "neos-3402294-bobin",
                                                  "neos-5104907-jarama",
                                                  "ns1644855",
                                                  "supportcase22",
                                                  "supportcase12",
                                                  "roi5alpha10n8",
                                                  "supportcase7",
                                                  "neos-4647030-tutaki",
                                                  "neos-5114902-kasavu",
                                                  "supportcase19",
                                                  "neos-5052403-cygnet",
                                                  "neos-2075418-temuka",
                                                  "neos-3402454-bohle",
                                                  "square41",
                                                  "square47"};

void print_result(const CSCMatrix<double>& matrix,
                  const std::vector<std::pair<size_t, size_t>>& result) {
  for (auto [i, j] : result) {
    if (i > j) {
      std::swap(i, j);
    }

    std::println("({}, {}):", i, j);

    auto xs = matrix.get_row(i);
    auto ys = matrix.get_row(j);

    for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      std::print("{:8}", x);
    }
    std::print("\n");

    for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      std::print("{:8}", y);
    }
    std::print("\n\n");
  }
}

void print_small_rows_cnt(const CSCMatrix<double>& matrix, size_t count) {
  auto [n, d] = matrix.shape();

  std::vector<size_t> counts(n, 0);
  for (size_t col = 0; col < d; ++col) {
    for (auto [row, value] : matrix.get_column(col)) {
      ++counts[row];
    }
  }

  std::vector<size_t> counts_counts(d + 1, 0);
  for (size_t row = 0; row < n; ++row) {
    ++counts_counts[counts[row]];
  }

  for (size_t i = 1; i <= count && i < counts_counts.size(); ++i) {
    std::println("  #{} - {}", i, counts_counts[i]);
  }
}

int main() {
  std::ofstream os("output.csv");
  std::println(os, "problem,found,considered");

  for (size_t i = 0; i < problems.size(); ++i) {
    const auto& problem_name = problems[i];

    auto matrix = get_problem_matrix(problem_name);

    std::println("{}/{} {} x {} ({})", i + 1, problems.size(),
                 matrix.shape().first, matrix.shape().second, problem_name);

    // print_small_rows_cnt(matrix, 6);

    auto seeker = seekers::Tables(4);
    auto result = seeker.seek(matrix);

    // if (result.size() != seekers::BruteForce(seekers::BruteForceHamming{2})
    //                          .seek(matrix)
    //                          .size()) {
    //   throw std::runtime_error("Wrong answer!");
    // }

    std::println(os, "{},{},{}", problem_name, result.size(),
                 seeker.get_stats().pairs_considered);

    //
    // std::println(output, "{},{},{}", problem_name, 4,
    //              seeker.get_stats().pairs_considered);
    //
    // output.flush();
    // }
  }

  // std::println("----------   tables   ----------");
  // print_result(matrix, result);
  //
  // std::println("---------- brute force ----------");
  // print_result(matrix, bf_result);

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
