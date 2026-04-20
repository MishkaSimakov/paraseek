#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "utils/Paths.h"

int main() {
  SCIP* scip = nullptr;

  // Create SCIP instance
  SCIP_CALL(SCIPcreate(&scip));

  // Include default plugins (presolve, heuristics, etc.)
  SCIP_CALL(SCIPincludeDefaultPlugins(scip));

  // Read MPS file
  auto path = paths::problem("problem1.mps");
  SCIP_CALL(SCIPreadProb(scip, path.c_str(), nullptr));

  // Transform problem (this triggers presolve)
  SCIP_CALL(SCIPtransformProb(scip));

  // Explicitly run presolve
  SCIP_CALL(SCIPpresolve(scip));

  // Print presolved problem statistics
  SCIP_CALL(SCIPprintStatistics(scip, nullptr));

  // If you want to stop here (only preprocessing)
  SCIPfree(&scip);

  return 0;
}
