/// Package config file. ONLY include in ".cc" files never in header files!

// std include
#include <string>

#ifndef ASKAP_ACCESSORS_H
#define ASKAP_ACCESSORS_H

  /// The name of the package
#define ASKAP_PACKAGE_NAME "synthesis"

/// askap namespace
namespace askap {
  /// @return version of the package
  std::string getAskapPackageVersion_synthesis();
}

  /// The version of the package
#define ASKAP_PACKAGE_VERSION askap::getAskapPackageVersion_synthesis()

#endif
