INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_MIMO_OFDM_JRC mimo_ofdm_jrc)

FIND_PATH(
    MIMO_OFDM_JRC_INCLUDE_DIRS
    NAMES mimo_ofdm_jrc/api.h
    HINTS $ENV{MIMO_OFDM_JRC_DIR}/include
        ${PC_MIMO_OFDM_JRC_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    MIMO_OFDM_JRC_LIBRARIES
    NAMES gnuradio-mimo_ofdm_jrc
    HINTS $ENV{MIMO_OFDM_JRC_DIR}/lib
        ${PC_MIMO_OFDM_JRC_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/mimo_ofdm_jrcTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MIMO_OFDM_JRC DEFAULT_MSG MIMO_OFDM_JRC_LIBRARIES MIMO_OFDM_JRC_INCLUDE_DIRS)
MARK_AS_ADVANCED(MIMO_OFDM_JRC_LIBRARIES MIMO_OFDM_JRC_INCLUDE_DIRS)
