#!/bin/bash
#===============================================================================
#         FILE: test/run_hook_interface_check.sh
#        USAGE: ./test/run_hook_interface_check.sh
#  DESCRIPTION: Prove that the per-stage user hook is CHECKED at compile time.
#
#               FUSS_solve and Explicit_Step declare their External_Function
#               dummy as procedure(external_function_if) -- a subroutine that
#               takes one type(FUSS_domain_type) argument. Before, the dummy
#               was `external`: any procedure was accepted, and hydra passed an
#               argument-less one for months without a diagnostic.
#
#               Two tiny programs are compiled against the preset build's
#               module directory (build/modules, where fuss_wrap_solve.mod
#               lives): one passes a procedure with the WRONG signature and
#               must fail to compile with an interface diagnostic; one passes
#               the right signature and must compile. A check that has never
#               been seen to fail is not a check, so the wrong one comes first.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FC=${FC:-gfortran}
MODDIR="$ROOT/build/modules"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/fuss-hook-check.XXXXXX")
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS hook interface check"
echo "  modules : $MODDIR"
echo "  work    : $WORK"
echo

if [ ! -f "$MODDIR/fuss_wrap_solve.mod" ]; then
  echo "  no $MODDIR/fuss_wrap_solve.mod -- build first: cmake --preset default && cmake --build build -j 8"
  exit 1
fi

cat > "$WORK/wrong.f90" <<'F90'
program hook_wrong
  use FUSS_Advanced_Types_m, only: FUSS_simulation_type
  use FUSS_Wrap_Solve,       only: FUSS_solve
  implicit none
  type(FUSS_simulation_type) :: simulation
  call FUSS_solve ( simulation, no_argument_hook )
contains
  subroutine no_argument_hook ()
  end subroutine no_argument_hook
end program hook_wrong
F90

cat > "$WORK/right.f90" <<'F90'
program hook_right
  use FUSS_Advanced_Types_m, only: FUSS_simulation_type, FUSS_domain_type
  use FUSS_Wrap_Solve,       only: FUSS_solve
  implicit none
  type(FUSS_simulation_type) :: simulation
  call FUSS_solve ( simulation, domain_hook )
contains
  subroutine domain_hook ( domain )
    type(FUSS_domain_type), intent(inout) :: domain
    if (.false.) domain%iter = domain%iter
  end subroutine domain_hook
end program hook_right
F90

# The wrong signature MUST be rejected, and rejected for the right reason: the
# diagnostic has to be about the procedure argument's interface, not some
# unrelated error. gfortran 12 words it "Interface mismatch in dummy procedure
# 'external_function'"; other compilers say "mismatch" too.
if "$FC" -c -I "$MODDIR" -o "$WORK/wrong.o" "$WORK/wrong.f90" > "$WORK/wrong.log" 2>&1; then
  bad "wrong signature COMPILED -- the hook is not checked (see $WORK/wrong.log)"
else
  if grep -qi "mismatch" "$WORK/wrong.log"; then
    pass "wrong signature rejected: $(grep -i -m1 'mismatch' "$WORK/wrong.log" | sed 's/^ *//')"
  else
    bad "wrong signature rejected, but not for an interface mismatch (see $WORK/wrong.log)"
  fi
fi

if "$FC" -c -I "$MODDIR" -o "$WORK/right.o" "$WORK/right.f90" > "$WORK/right.log" 2>&1; then
  pass "conforming signature compiles"
else
  bad "conforming signature FAILED to compile (see $WORK/right.log)"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "hook interface check passed"
  rm -rf "$WORK"
else
  echo "$fail hook interface check(s) FAILED -- see $WORK"
  exit 1
fi
