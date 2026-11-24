# Each invocation owns a scratch directory. Exact status 1 rejects both signals
# and accidental usage success; stderr must name the inaccessible input. Pitfall:
# WILL_FAIL alone also accepts crashes and cannot establish this error contract.
file(MAKE_DIRECTORY "${SCRATCH}")
file(GLOB before "${SCRATCH}/*")
if(before)
    message(FATAL_ERROR "Expected an empty scratch directory: ${SCRATCH}")
endif()
set(input "${SCRATCH}/missing.${EXTENSION}")
set(args "${input}" 10)
if(TIME_DRIVER)
    list(APPEND args 0)
endif()
execute_process(COMMAND "${DRIVER}" ${args}
    WORKING_DIRECTORY "${SCRATCH}"
    RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr
    TIMEOUT 3)
if(NOT "${status}" STREQUAL "1")
    message(FATAL_ERROR "Expected status 1, got ${status}: ${stderr}")
endif()
string(FIND "${stderr}" "${input}" filename_pos)
string(FIND "${stderr}" "Cannot read moment file:" message_pos)
if(filename_pos EQUAL -1 OR message_pos EQUAL -1)
    message(FATAL_ERROR "Missing filename or diagnostic on stderr: ${stderr}")
endif()
file(GLOB after "${SCRATCH}/*")
if(after)
    message(FATAL_ERROR "Driver wrote files on input failure: ${after}")
endif()
file(REMOVE_RECURSE "${SCRATCH}")
