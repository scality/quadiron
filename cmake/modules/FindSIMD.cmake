#
# Check if SSE4.1 or AVX2 instructions are available on the machine where the
# project is compiled.
# These commands are extracted and modified from
# https://github.com/VcDevel/Vc/blob/master/cmake/FindSSE.cmake
#
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
   EXEC_PROGRAM(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)

   STRING(REGEX REPLACE "^.*(sse4_1).*$" "\\1" SSE_THERE ${CPUINFO})
   STRING(COMPARE EQUAL "sse4_1" "${SSE_THERE}" SSE41_TRUE)
   if (SSE41_TRUE)
      set(SSE4_1_FOUND true)
   else (SSE41_TRUE)
      set(SSE4_1_FOUND false)
   endif (SSE41_TRUE)

   STRING(REGEX REPLACE "^.*(avx2).*$" "\\1" SSE_THERE ${CPUINFO})
   STRING(COMPARE EQUAL "avx2" "${SSE_THERE}" AVX2_TRUE)
   if (AVX2_TRUE)
      set(AVX2_FOUND true)
   else (AVX2_TRUE)
      set(AVX2_FOUND false)
   endif (AVX2_TRUE)
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
   EXEC_PROGRAM("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE
      CPUINFO)

   STRING(REGEX REPLACE "^.*(SSE4.1).*$" "\\1" SSE_THERE ${CPUINFO})
   STRING(COMPARE EQUAL "SSE4.1" "${SSE_THERE}" SSE41_TRUE)
   if (SSE41_TRUE)
      set(SSE4_1_FOUND true)
   else (SSE41_TRUE)
      set(SSE4_1_FOUND false)
   endif (SSE41_TRUE)

   STRING(REGEX REPLACE "^.*(AVX1.0).*$" "\\1" SSE_THERE ${CPUINFO})
   STRING(COMPARE EQUAL "AVX1.0" "${SSE_THERE}" AVX2_TRUE)
   if (AVX2_TRUE)
      set(AVX2_FOUND true)
   else (AVX2_TRUE)
      set(AVX2_FOUND false)
   endif (AVX2_TRUE)
else()
   set(SSE4_1_FOUND true)
   set(AVX2_FOUND false)
endif()

if(NOT SSE4_1_FOUND)
  MESSAGE(STATUS "Could not find hardware support for SSE4.1 on this machine.")
endif(NOT SSE4_1_FOUND)

if(NOT AVX2_FOUND)
  MESSAGE(STATUS "Could not find hardware support for AVX2 on this machine.")
endif(NOT AVX2_FOUND)

#
# Determine highest level for supported SIMD
#
if (AVX2_FOUND)
  set(HIGHEST_SIMD "AVX2")
  set(SIMD_LIST "OFF" "SSE4" "AVX2")
elseif (SSE4_1_FOUND)
  set(HIGHEST_SIMD "SSE4")
  set(SIMD_LIST "OFF" "SSE4")
else()
  set(HIGHEST_SIMD "OFF")
  set(SIMD_LIST "OFF")
endif()

#
# Determine suitable SIMD level
#
if(NTTEC_USE_SIMD STREQUAL "AVX2")
  if(NOT AVX2_FOUND)
    set(NTTEC_USE_SIMD ${HIGHEST_SIMD})
    MESSAGE(WARNING "AVX2 enabled but not supported: fallback on ${NTTEC_USE_SIMD}.")
  endif()
elseif(NTTEC_USE_SIMD STREQUAL "SSE4")
  if(NOT SSE4_1_FOUND)
    set(NTTEC_USE_SIMD ${HIGHEST_SIMD})
    MESSAGE(WARNING "SSE4 enabled but not supported: fallback on ${NTTEC_USE_SIMD}.")
  endif()
elseif(NTTEC_USE_SIMD STREQUAL "ON")
  # if level is not specified, choose highest value
  set(NTTEC_USE_SIMD ${HIGHEST_SIMD})
elseif(NOT NTTEC_USE_SIMD)
  set(NTTEC_USE_SIMD "OFF")
else()
  MESSAGE(WARNING "Given ${NTTEC_USE_SIMD} is not recognised for SIMD option. Set it OFF.")
  set(NTTEC_USE_SIMD "OFF")
endif()

#
# Add definition to distinguish SSE and AVX
#
if(NTTEC_USE_SIMD)
  add_definitions(-DNTTEC_USE_${NTTEC_USE_SIMD})
endif()
