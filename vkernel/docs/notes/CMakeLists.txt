# The studies of the notes of vkernel/docs/notes: one target each, not
# part of all, writing their outputs (data, tables, figures) under
# <build>/vkernel/soa/notes/<study>. The figures of the documents are
# copies of these outputs.
#
#   cmake --build <build> --target vkernel_study_polytope_contact_manifolds
#   cmake --build <build> --target vkernel_study_time_integrators
#
# The figures need matplotlib, which is in the development environment of
# guix.scm: without it the studies still run, with the data and the tables
# only.

set(_notes_src ${CMAKE_CURRENT_SOURCE_DIR})
set(_notes_bin ${CMAKE_CURRENT_BINARY_DIR})

# polytope-contact-manifolds.md: vnative scenes through the python runner
if(TARGET _nonos AND WITH_PYB11_WRAPPER)
  set(_libs)
  foreach(_component IN LISTS COMPONENTS)
    if(TARGET ${_component})
      list(APPEND _libs $<TARGET_FILE_DIR:${_component}>)
    endif()
  endforeach()
  foreach(_target IN ITEMS _nonos CompactNSearch cuNSearch)
    if(TARGET ${_target})
      list(APPEND _libs $<TARGET_FILE_DIR:${_target}>)
    endif()
  endforeach()
  string(JOIN ":" _libs ${_libs})
  add_custom_target(vkernel_study_polytope_contact_manifolds
    COMMAND ${Python_EXECUTABLE}
            ${_notes_src}/polytope-contact-manifolds/scripts/run.py
            --build-dir ${CMAKE_BINARY_DIR}
            --out ${_notes_bin}/polytope-contact-manifolds
            --nonos-module $<TARGET_FILE:_nonos>
            --library-path "${_libs}"
    DEPENDS _nonos
    COMMENT "study: polytope contact manifolds"
    USES_TERMINAL
    VERBATIM)
endif()

# time-integrators-mj-at-dem.md: nscd_comparisons, then the tables and plots
if(TARGET nscd_comparisons)
  set(_out ${_notes_bin}/time-integrators-mj-at-dem)
  add_custom_target(vkernel_study_time_integrators
    COMMAND ${CMAKE_COMMAND} -E make_directory ${_out}
    COMMAND $<TARGET_FILE:nscd_comparisons> ${_out}
    COMMAND ${Python_EXECUTABLE}
            ${_notes_src}/time-integrators-mj-at-dem/scripts/nscd_comparisons_plot.py
            ${_out}
    DEPENDS nscd_comparisons
    COMMENT "study: Moreau-Jean, Anitescu-Tasora, DEM"
    USES_TERMINAL
    VERBATIM)
endif()
