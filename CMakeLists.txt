cmake_minimum_required (VERSION 3.8)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("Sp3Resample")

add_executable (Sp3Resample "main.c"  "sp3.h" "sp3.c"   "interp.h" "interp.c" "linsolve.h" "linsolve.c" "astrodyn.h" "astrodyn.c" "error.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Sp3Resample PROPERTY CXX_STANDARD 20)
endif()
