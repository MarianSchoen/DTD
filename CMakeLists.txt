cmake_minimum_required (VERSION 3.0)
set(PROJECT "DTD_CXX")
project (${PROJECT} VERSION "0.1.0")

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
    set (CMAKE_CXX_STANDARD_REQUIRED ON)
  endif ()
endmacro(use_cxx11)

use_cxx11()

add_library( dtd STATIC src/utils.cc src/models.cc src/fista.hpp src/utils.hpp src/models.hpp) #add header files here
set_target_properties(dtd PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

add_library( dtd_jl SHARED src/interface_jl.cc src/interface_jl.hpp )
target_link_libraries(dtd_jl dtd)
set_target_properties(dtd_jl PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Werror -Wshadow -Wnon-virtual-dtor -pedantic")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror -Wshadow -Wnon-virtual-dtor -pedantic")
SET(CMAKE_CXX_FLAGS_PROF "-O3 -pg -g -Wall -Werror -Wshadow -Wnon-virtual-dtor -pedantic")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(dtd Eigen3::Eigen)

#??? 
find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(PUBLIC my_app "-pthread")
endif()

install(TARGETS dtd
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
