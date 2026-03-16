#include <benchmark/benchmark.h>

#include "ProblemMatrix.h"
#include "seekers/BruteForce.h"

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

static void benchmark_brute_force(benchmark::State& state) {
  const std::string& problem_name = problems[state.range(0)];

  state.PauseTiming();
  auto matrix = get_problem_matrix(problem_name);
  state.ResumeTiming();

  for (auto _ : state) {
    seekers::BruteForce(3).seek(matrix);
  }

  state.SetComplexityN(matrix.nonzero_count());
}

int main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);

  benchmark::RegisterBenchmark("BruteForce", benchmark_brute_force)
      ->DenseRange(0, problems.size() - 1)
      ->Complexity();

  benchmark::RunSpecifiedBenchmarks();
}
