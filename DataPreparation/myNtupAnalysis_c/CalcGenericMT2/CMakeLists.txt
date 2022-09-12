################################################################################
# Package: CalcGenericMT2
################################################################################

# Declare the package name:
atlas_subdir( CalcGenericMT2 )

# atlas_depends_on_subdirs is deprecated in release22
if($ENV{AtlasVersion} MATCHES "^21")
  # Declare the package's dependencies:
  atlas_depends_on_subdirs( PUBLIC
                          Event/xAOD/xAODBase
                          Event/xAOD/xAODMissingET
)
endif()

# Find the needed external(s):
find_package( ROOT COMPONENTS Core RIO Hist Tree )

# build a dictionary for the library
atlas_add_root_dictionary ( CalcGenericMT2Lib CalcGenericMT2DictSource
                            ROOT_HEADERS CalcGenericMT2/MT2_*.h Root/LinkDef.h
                            EXTERNAL_PACKAGES ROOT
)

# build a shared library
atlas_add_library( CalcGenericMT2Lib CalcGenericMT2/*.h ${CalcGenericMT2DictSource}
                   PUBLIC_HEADERS CalcGenericMT2
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
                   xAODBase xAODMissingET
)
