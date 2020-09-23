###################################################################################
#
# Vessel tree construction
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

if(enable_GLOG)
  set(compile_options_tmp ${base_compile_options} " -D_GLOG")
elseif()
  set(compile_options_tmp ${base_compile_options})
endif()

string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

set(main_target "minimumComplianceProblem")

add_executable(${main_target} ${main_src})
target_sources(${main_target}
  PRIVATE
    main_pardiso.cpp
)
target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

target_link_libraries(${main_target} PRIVATE
  MinimumCompliance
  FEM
  LS
  Base
  ${GLOG_LIBRARY}
  ${TP_LIBRARY}
  ${MKL_LIBRARIES}
)

target_include_directories(${main_target} PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib/FEM
  ${PROJECT_SOURCE_DIR}/src/shapeAndTopologyOptimization/shapeOptimization/minimumComplianceProblem
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})