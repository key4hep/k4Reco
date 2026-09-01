# - Locate FastJet library (including fastjet-contrib)
# Defines:
#
#  FASTJET_FOUND
#  FASTJET_INCLUDE_DIR
#  FASTJET_INCLUDE_DIRS (not cached)
#  FASTJET_LIBRARY
#  FASTJET_LIBRARIES (not cached)
#  FASTJET_LIBRARY_DIRS (not cached)

find_path(FASTJET_INCLUDE_DIR fastjet/version.hh
          HINTS $ENV{FASTJET_ROOT_DIR}/include ${FASTJET_ROOT_DIR}/include)

# fastjet-contrib provides the ValenciaPlugin. It is a separate package that may be
# installed either into the FastJet prefix or into one of its own
find_path(FASTJET_CONTRIB_INCLUDE_DIR fastjet/contrib/ValenciaPlugin.hh
          HINTS $ENV{FJCONTRIB_ROOT_DIR}/include ${FJCONTRIB_ROOT_DIR}/include
                $ENV{FASTJET_ROOT_DIR}/include ${FASTJET_ROOT_DIR}/include
                ${FASTJET_INCLUDE_DIR})

find_library(FASTJET_LIBRARY NAMES fastjet
             HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

find_library(FASTJETPLUGINS_LIBRARY NAMES fastjetplugins
             HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

# Needed by fastjet-contrib, which uses the background estimators and Recluster
find_library(FASTJETTOOLS_LIBRARY NAMES fastjettools
             HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

find_library(FASTJETCONTRIB_LIBRARY NAMES fastjetcontribfragile fastjetcontrib
             HINTS $ENV{FJCONTRIB_ROOT_DIR}/lib ${FJCONTRIB_ROOT_DIR}/lib
                   $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

find_library(SISCONE_LIBRARY NAMES siscone
             HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

find_library(SISCONE_SPHERICAL_LIBRARY NAMES siscone_spherical
             HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

# handle the QUIETLY and REQUIRED arguments and set FASTJET_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FastJet DEFAULT_MSG FASTJET_INCLUDE_DIR FASTJET_LIBRARY
                                  FASTJET_CONTRIB_INCLUDE_DIR FASTJETCONTRIB_LIBRARY
                                  FASTJETTOOLS_LIBRARY)

mark_as_advanced(FASTJET_FOUND FASTJET_INCLUDE_DIR FASTJET_LIBRARY FASTJET_CONTRIB_INCLUDE_DIR
                 FASTJETPLUGINS_LIBRARY FASTJETTOOLS_LIBRARY FASTJETCONTRIB_LIBRARY
                 SISCONE_LIBRARY SISCONE_SPHERICAL_LIBRARY)

set(FASTJET_INCLUDE_DIRS ${FASTJET_INCLUDE_DIR} ${FASTJET_CONTRIB_INCLUDE_DIR})
list(REMOVE_DUPLICATES FASTJET_INCLUDE_DIRS)
# Order matters: the contrib library has unresolved references into fastjettools, and with
# --as-needed the linker only keeps a library that resolves symbols still undefined when it
# is reached. Dependents must therefore come before their dependencies
set(FASTJET_LIBRARIES ${FASTJETCONTRIB_LIBRARY} ${FASTJETTOOLS_LIBRARY} ${FASTJETPLUGINS_LIBRARY}
                      ${SISCONE_SPHERICAL_LIBRARY} ${SISCONE_LIBRARY} ${FASTJET_LIBRARY})
get_filename_component(FASTJET_LIBRARY_DIRS ${FASTJET_LIBRARY} PATH)
