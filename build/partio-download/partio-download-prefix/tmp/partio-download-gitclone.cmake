
if(NOT "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitinfo.txt" IS_NEWER_THAN "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone  "https://github.com/wdas/partio.git" "partio-src"
    WORKING_DIRECTORY "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/wdas/partio.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout v1.1.0 --
  WORKING_DIRECTORY "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'v1.1.0'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  submodule update --recursive --init 
  WORKING_DIRECTORY "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitinfo.txt"
    "/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/partio-download/partio-download-prefix/src/partio-download-stamp/partio-download-gitclone-lastrun.txt'")
endif()

