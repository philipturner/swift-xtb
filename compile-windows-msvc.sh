# In this attempt, try compiling with Intel oneAPI compilers instead of GCC.
# While at the same time, integrating with existing system executables instead
# of a containerized environment. Avoid choices that will cause symbol
# conflicts between two libraries trying to act as MSVC.
#
# Start by compiling a single Fortran file with ifx/icx.
#
# If this doesn't work out, the last resort is MSYS2 + runtime symbol loading.
