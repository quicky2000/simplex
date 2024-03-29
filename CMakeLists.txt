cmake_minimum_required(VERSION 3.8)
project(simplex)

set(CMAKE_CXX_STANDARD 11)

# List source files of project
set(MY_SOURCE_FILES
    include/simplex_array.h
    include/simplex_array_base.h
    include/simplex_listener.h
    include/simplex_map.h
    include/simplex_solver_base.h
    include/simplex_solver_integer_base.h
    include/simplex_solver_integer.h
    include/simplex_solver_integer_ppcm.h
    include/simplex_solver.h
    include/simplex_identity_solver.h
    include/my_matrix.h
    include/my_square_matrix.h
    include/equation_system.h
    src/test_my_matrix.cpp
    src/test_my_square_matrix.cpp
    src/test_equation_system.cpp
    include/simplex_listener_target_if.h
    include/simplex_solver_glpk.h
    include/simplex_listener_if.h
    src/simplex_solver_glpk.cpp
   )

# List header directories in project
set(MY_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

#Create dependancy list
set(DEPENDANCY_LIST "quicky_utils")
LIST(APPEND DEPENDANCY_LIST "quicky_exception")
LIST(APPEND DEPENDANCY_LIST "EXT_glpk")

#------------------------------
#- Generic part
#------------------------------
set(DEPENDANCY_OBJECTS "")
set(LINKED_LIBRARIES "")
#iterate on dependancies
foreach(DEPENDANCY_ITEM IN ITEMS ${DEPENDANCY_LIST})
    message("${PROJECT_NAME} need dependancy ${DEPENDANCY_ITEM}")
    set(IS_IMPORTED FALSE)
    #Get information from dependancies if not already done
    if(NOT TARGET ${DEPENDANCY_ITEM})
        set(SUB_DEPENDANCY_OBJECTS "")
        set(SUB_LINKED_LIBRARIES "")
        message("${PROJECT_NAME}-------------------before |${DEPENDANCY_OBJECTS}|")
        message("${PROJECT_NAME} include dependancy ${DEPENDANCY_ITEM}")
        add_subdirectory($ENV{QUICKY_REPOSITORY}/${DEPENDANCY_ITEM} ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDANCY_ITEM})
        message("${PROJECT_NAME}-------------------after_sub |${SUB_DEPENDANCY_OBJECTS}|")
        get_target_property(IS_IMPORTED ${DEPENDANCY_ITEM} IMPORTED)
        if(NOT ${IS_IMPORTED})
            list(APPEND DEPENDANCY_OBJECTS "$<TARGET_OBJECTS:${DEPENDANCY_ITEM}>" ${SUB_DEPENDANCY_OBJECTS})
        else()
            list(APPEND LINKED_LIBRARIES ${DEPENDANCY_ITEM})
        endif()
        list(APPEND LINKED_LIBRARIES ${SUB_LINKED_LIBRARIES})
        message("${PROJECT_NAME}-------------------after |${DEPENDANCY_OBJECTS}|")
    endif()

    # Complete include directories with directories coming from dependancy because of transitivity issues with Object libraries
    if(NOT ${IS_IMPORTED})
        set(MY_INCLUDE_DIRECTORIES ${MY_INCLUDE_DIRECTORIES} $<TARGET_PROPERTY:${DEPENDANCY_ITEM},INCLUDE_DIRECTORIES>)
    else(NOT ${IS_IMPORTED})
        set(MY_INCLUDE_DIRECTORIES ${MY_INCLUDE_DIRECTORIES} $<TARGET_PROPERTY:${DEPENDANCY_ITEM},INTERFACE_INCLUDE_DIRECTORIES>)
    endif()
endforeach(DEPENDANCY_ITEM)


#Prepare targets
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(IS_DIRECTORY ${HAS_PARENT})
    message("Declare library ${PROJECT_NAME}")
    add_library(${PROJECT_NAME} OBJECT ${MY_SOURCE_FILES})
    set(SUB_DEPENDANCY_OBJECTS ${DEPENDANCY_OBJECTS} PARENT_SCOPE)
    set(SUB_LINKED_LIBRARIES ${LINKED_LIBRARIES} PARENT_SCOPE)
else()
    #set(CMAKE_VERBOSE_MAKEFILE ON)
    add_executable(${PROJECT_NAME} ${MY_SOURCE_FILES} ${DEPENDANCY_OBJECTS} src/main_${PROJECT_NAME}.cpp)
    target_compile_definitions(${PROJECT_NAME} PRIVATE -D$<UPPER_CASE:${PROJECT_NAME}>_SELF_TEST)
    target_compile_options(${PROJECT_NAME} PUBLIC -Wall -pedantic -g -O0)
    target_link_libraries(${PROJECT_NAME} ${LINKED_LIBRARIES})
endif()

target_include_directories(${PROJECT_NAME} PUBLIC ${MY_INCLUDE_DIRECTORIES})

foreach(DEPENDANCY_ITEM IN ITEMS ${DEPENDANCY_LIST})
    add_dependencies(${PROJECT_NAME} ${DEPENDANCY_ITEM})
endforeach(DEPENDANCY_ITEM)

#Force use of -std=c++11 instead of -std=gnu++11
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

#EOF
