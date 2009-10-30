#include "ControlPlotsComparison.h"

#include <string>
#include <vector>

int main()
{
  std::vector<std::string> files;
  std::vector<std::string> descr;

  // Add files with controlplots to be compared
  files.push_back("test/controlplots_1000.root");
  files.push_back("test/controlplots_-1.root");

  // Add descriptions specifying the different
  // files; the descriptions will appear in the
  // legend of the comparison plots
  descr.push_back("1000 events");
  descr.push_back("All events");

  // Create comparison plots
  TControlPlotsComparison comp;
  comp.CompareControlPlots(files,descr);

  return 0;
}
