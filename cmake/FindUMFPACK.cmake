IF(Umf_INCLUDE_DIRS)

  FIND_PATH(UMF_INCLUDE_DIR umfpack.h ${Umf_INCLUDE_DIRS})
  FIND_LIBRARY(UMF_LIBRARY umfpack ${Umf_LIBRARY_DIRS})

ELSE(Umf_INCLUDE_DIRS)

  SET(TRIAL_LIBRARY_PATHS
    /usr/lib 
    /usr/local/lib
    /opt/lib
    /sw/lib
    )

  SET(TRIAL_INCLUDE_PATHS
    /usr/include
    /usr/include/suitesparse
    /opt/include
    /usr/local/include
    /sw/include
    )

  FIND_LIBRARY(UMF_LIBRARY umfpack ${TRIAL_LIBRARY_PATHS})
  FIND_PATH(UMF_INCLUDE_DIR umfpack.h ${TRIAL_INCLUDE_PATHS} )

ENDIF(Umf_INCLUDE_DIRS)

IF(UMF_INCLUDE_DIR AND UMF_LIBRARY)
  SET(FOUND_UMF 1 CACHE BOOL "Found umfpack library")
  SET(UMF_LIBRARY ${UMF_LIBRARY} -lumfpack)
ELSE(UMF_INCLUDE_DIR AND UMF_LIBRARY)
  SET(FOUND_UMF 0 CACHE BOOL "Not fount umf library")
ENDIF(UMF_INCLUDE_DIR AND UMF_LIBRARY)

MARK_AS_ADVANCED(
  UMF_INCLUDE_DIR 
  UMF_LIBRARY 
  FOUND_UMF
  )
