# - Try to find the TALPA (ELTOPO) includes
# Once done this will define
#
#  TALPA_FOUND - system has ELTOPO/COMMON
#  LIBTALPA_INCLUDE_DIR - the ELTOPO/TALPA include directory
if(NOT LIBTALPA_FOUND)

FIND_PATH(LIBTALPA_INCLUDE_DIR iomesh.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_SOURCE_DIR}/../eltopo/talpa
   /usr/include
   /usr/local/include
)

# Includes must be found
if(NOT LIBTALPA_INCLUDE_DIR)
  set(LIBTALPA_INCLUDE_DIR FALSE)
  message(FATAL_ERROR "could NOT find LIBTALPA_INCLUDE_DIR")
endif(NOT LIBTALPA_INCLUDE_DIR)

if(LIBTALPA_INCLUDE_DIR)
   set(LIBTALPA_FOUND TRUE)
   set(LIBTALPA_INCLUDE_DIRS ${LIBTALPA_INCLUDE_DIR})
endif(LIBTALPA_INCLUDE_DIR)

endif(NOT LIBTALPA_FOUND)