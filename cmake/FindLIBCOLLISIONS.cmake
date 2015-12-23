# - Try to find the LIBCOLLISIONS library
# Once done this will define
#
#  COLLISIONS_FOUND - system has LIBCOLLISIONS
#  LIBCOLLISIONS_INCLUDE_DIR - the LIBCOLLISIONS include directory
#  LIBCOLLISIONS_LYBRARY - the LIBCOLLISIONS library
if(NOT LIBCOLLISIONS_FOUND)

FIND_PATH(LIBCOLLISIONS_INCLUDE_DIR CTCD.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_SOURCE_DIR}/../collisiondetection/include/
   ${PROJECT_SOURCE_DIR}/../../collisiondetection/include/
   ${PROJECT_SOURCE_DIR}/../../../collisiondetection/include/
   ${PROJECT_SOURCE_DIR}/../code/collisiondetection/include/
   ${PROJECT_SOURCE_DIR}/../../code/collisiondetection/include/
   ${PROJECT_SOURCE_DIR}/../../../code/collisiondetection/include/
   /usr/include
   /usr/local/include
)

FIND_PATH(LIBCOLLISIONS_MAIN_DIR Makefile
   ${PROJECT_SOURCE_DIR}/../collisiondetection/
   ${PROJECT_SOURCE_DIR}/../../collisiondetection/
   ${PROJECT_SOURCE_DIR}/../../../collisiondetection/
   ${PROJECT_SOURCE_DIR}/../code/collisiondetection/
   ${PROJECT_SOURCE_DIR}/../../code/collisiondetection/
   ${PROJECT_SOURCE_DIR}/../../../code/collisiondetection/
)

# Includes must be found
if(NOT LIBCOLLISIONS_INCLUDE_DIR)
  set(LIBCOLLISIONS_INCLUDE_DIR FALSE)
  message(FATAL_ERROR "could NOT find LIBCOLLISIONS_INCLUDE_DIR")
endif(NOT LIBCOLLISIONS_INCLUDE_DIR)

if(LIBCOLLISIONS_INCLUDE_DIR)
   set(LIBCOLLISIONS_FOUND TRUE)
   set(LIBCOLLISIONS_INCLUDE_DIRS ${LIBCOLLISIONS_INCLUDE_DIR})
endif(LIBCOLLISIONS_INCLUDE_DIR)

FIND_LIBRARY(LIBCOLLISIONS_LIBRARY NAME collisions PATHS ${LIBCOLLISIONS_MAIN_DIR}/bin)

# Main library must be found
if(NOT LIBCOLLISIONS_LIBRARY)
  set(LIBCOLLISIONS_FOUND FALSE)
  message(FATAL_ERROR "could NOT find LIBCOLLISIONS")
endif(NOT LIBCOLLISIONS_LIBRARY)




if(LIBCOLLISIONS_FOUND)
   if(NOT LIBCOLLISIONS_FIND_QUIETLY)
      message(STATUS "Found LIBCOLLISIONS: ${LIBCOLLISIONS_LIBRARY}")
   endif(NOT LIBCOLLISIONS_FIND_QUIETLY)
else(LIBCOLLISIONS_FOUND)
   if(LIBCOLLISIONS_FIND_REQUIRED)
      message(FATAL_ERROR "could NOT find LIBCOLLISIONS")
   endif(LIBCOLLISIONS_FIND_REQUIRED)
endif(LIBCOLLISIONS_FOUND)

endif(NOT LIBCOLLISIONS_FOUND)