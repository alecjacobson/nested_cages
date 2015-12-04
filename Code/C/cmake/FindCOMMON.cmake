# - Try to find the COMMON (ELTOPO) includes
# Once done this will define
#
#  COMMON_FOUND - system has ELTOPO/COMMON
#  LIBCOMMON_INCLUDE_DIR - the ELTOPO/COMMON include directory
if(NOT LIBCOMMON_FOUND)

FIND_PATH(LIBCOMMON_INCLUDE_DIR commonoptions.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_SOURCE_DIR}/../eltopo/common
   /usr/include
   /usr/local/include
)

# Includes must be found
if(NOT LIBCOMMON_INCLUDE_DIR)
  set(LIBCOMMON_INCLUDE_DIR FALSE)
  message(FATAL_ERROR "could NOT find LIBCOMMON_INCLUDE_DIR")
endif(NOT LIBCOMMON_INCLUDE_DIR)

if(LIBCOMMON_INCLUDE_DIR)
   set(LIBCOMMON_FOUND TRUE)
   set(LIBCOMMON_INCLUDE_DIRS ${LIBCOMMON_INCLUDE_DIR})
endif(LIBCOMMON_INCLUDE_DIR)

endif(NOT LIBCOMMON_FOUND)