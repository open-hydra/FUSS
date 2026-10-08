#!/bin/bash -
#===============================================================================
#
#          FILE: FUSS.sh
#
#         USAGE: run "./FUSS.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run FUSS
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     compile accordingly with the preset file"
  echo " solve                 -->     run FUSS"
  echo " kill                  -->     kill the process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -m | --mpi <n>        -->     launch solver with <n> MPI ranks"
  echo " -p | --parallel <n>   -->     launch solver with <n> OMP threads per MPI rank"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../..
MASTER=$MASTERDIR/bin/FUSS
LOCAL=./bin/FUSS

# Default Options
BG=0
NMPI=1
NTHREADS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    -b | --background)
        #echo " -> Background"
        BG=1
        shift
        ;;
    -m | --mpi)
        NMPI=$2
        shift 2
        ;;
    -p | --parallel)
        #echo " -> OpenMP($2)"
        NTHREADS=$2
        shift 2
        ;;
    -h | --help)
        print_usage
        shift
        ;;
    * )
      break
      ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi
   
if [[ $1 == solve ]]; then
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  # Run the solver
  if [[ $NMPI == 1 ]]; then
    if [[ $BG == 0 ]]; then
      $LOCAL
      RC=$?
    else
      $LOCAL 2>errors_file >logfile &
      echo $! > .ID
    fi
  else
    if [[ $BG == 0 ]]; then
      mpirun -np $NMPI --map-by socket --bind-to socket $LOCAL
      RC=$?
    else
      mpirun -np $NMPI --map-by socket --bind-to socket $LOCAL 2>errors_file >logfile &
      echo $! > .ID
    fi
  fi
fi

if [[ $1 == kill ]]; then
read PID < .ID && kill $PID
fi

# Propagate the solver's exit status.
#
# Without this the script always exits 0, because the last thing it evaluates is
# the `kill` test above. That silently defeated every abort code the solver
# produces: test/run_regression.sh branches on exit status and recorded a case
# that had refused to run as 'ok', fingerprinting whatever partial output it
# managed to write. Codes: 1 input validation, 2 geometry/mesh, 3 numerical.
#
# Background runs cannot report a status here -- the solver is still running --
# so RC stays 0 for them, by design.
exit ${RC:-0}
