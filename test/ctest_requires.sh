#!/bin/bash
#===============================================================================
#         FILE: test/ctest_requires.sh
#        USAGE: ./test/ctest_requires.sh REQUIREMENT... -- COMMAND [ARG...]
#  DESCRIPTION: Runs COMMAND (exec) when every REQUIREMENT holds; otherwise
#               prints one SKIP line with the reason and exits 77, the
#               SKIP_RETURN_CODE of the gate tests of test/CMakeLists.txt, so
#               ctest reports the test Not Run instead of Failed. By hand the
#               gate scripts are run directly and keep their own refusals.
#
#               REQUIREMENT
#                 impi              I_MPI_ROOT is set and $I_MPI_ROOT/bin/mpirun is
#                                   executable (the Intel MPI the MPI gates use)
#                 file=PATH         PATH exists
#                 fc=PATH[:VERSION] the compiler PATH is executable and, with
#                                   VERSION, `PATH -dumpfullversion` prints it:
#                                   gfortran reads only the .mod files its own
#                                   version wrote
#===============================================================================
set -u

skip () { echo "SKIP: $1"; exit 77; }

while [ $# -gt 0 ] && [ "$1" != "--" ]; do
  case $1 in
    impi)
      [ -n "${I_MPI_ROOT:-}" ] || skip "I_MPI_ROOT is not set (the Intel MPI of the mpi preset)"
      [ -x "$I_MPI_ROOT/bin/mpirun" ] || skip "no executable $I_MPI_ROOT/bin/mpirun" ;;
    file=*)
      [ -e "${1#file=}" ] || skip "no ${1#file=}" ;;
    fc=*)
      spec=${1#fc=}; fc=${spec%%:*}; want=""
      [ "$spec" != "$fc" ] && want=${spec#*:}
      [ -x "$fc" ] || command -v "$fc" > /dev/null 2>&1 || skip "no Fortran compiler $fc"
      if [ -n "$want" ]; then
        have=$("$fc" -dumpfullversion 2> /dev/null)
        [ "$have" = "$want" ] || skip "$fc is version ${have:-unknown}, the build tree was configured with $want"
      fi ;;
    *) echo "ctest_requires.sh: unknown requirement $1" >&2; exit 2 ;;
  esac
  shift
done
[ "${1:-}" = "--" ] || { echo "usage: $0 REQUIREMENT... -- COMMAND [ARG...]" >&2; exit 2; }
shift
[ $# -gt 0 ] || { echo "usage: $0 REQUIREMENT... -- COMMAND [ARG...]" >&2; exit 2; }
exec "$@"
