#include "scip/scip.h"

#include "scip/scipdefplugins.h"
#include "utils/Paths.h"

int main(int argc, char** argv) {
  SCIP* scip = NULL;

  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);

  // Read MPS
  auto path = paths::problem("square47.mps");
  SCIPreadProb(scip, path.c_str(), NULL);

  // Presolve
  SCIPpresolve(scip);

  // Write transformed problem
  auto output_path = paths::problem("presolved_square47.mps");
  SCIPwriteTransProblem(scip, output_path.c_str(), "mps", FALSE);

  // Cleanup
  SCIPfree(&scip);

  return 0;
}
