#!/bin/sh
# ChainHash-256 certificate suite, on the as-built reference (test/256/ph_ref.h) at every block size:
#   field: Rabin irreducibility of x^256+x^10+x^5+x^2+1 (with negative controls), and the C field product
#          against an independent Python product;
#   level 1 (cert_pk.c, power key): the as-built block value, as a polynomial in the key s over GF(2^256), equals
#          the exponent-class expansion sum x_i y_i + sum (x_i s^2i + y_i s^(2i-1)) + sum s^(4i-1) (interpolated,
#          then checked at extra nodes); adversarial differences give nonzero difference polynomials of degree
#          <= 2 PH_M, PH_M = 8 (quick) and 64 (the released block size); a colliding-exponent derivation is rejected;
#   level 2 (cert_2l.c): the region polynomial equals Lemma A's expansion for q = 1..8, the outer polynomial
#          in z, the <= 1-block path; a negative control (colliding exponents) must be rejected;
#   bound: ph256_bounds.py, N_2L(L) <= 2L (equality only at L = 1; score 255).
# usage: test/256/certs.sh   (CC=... to choose the compiler). Exit status 0 iff everything passes.
set -u
here=$(cd "$(dirname "$0")" && pwd); CC=${CC:-cc}; fail=0
tmp=$(mktemp -d "${TMPDIR:-/tmp}/ch256cert.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
case $(uname -m) in x86_64|amd64) HW="-mpclmul -msse4.1 -DCH256_HW_CLMUL";; *) HW="-DCH256_HW_CLMUL";; esac
build() { $CC -O2 $HW -I"$here" "$@"; }
# check <lines to show> <command...>: run, show the last lines of its output, fail on a nonzero exit status
check() { n=$1; shift; "$@" > "$tmp/out" 2>&1; rc=$?; tail -n "$n" "$tmp/out"; [ $rc = 0 ] || { echo "FAIL (exit $rc): $*"; fail=1; }; }
# verdict <pattern>: the last check's output must contain the pattern and no FAIL / ACCEPTED line
verdict() { grep -q "$1" "$tmp/out" && ! grep -qE 'FAIL|ACCEPTED' "$tmp/out" || { echo "FAIL: no '$1' verdict"; fail=1; }; }
check 3 python3 "$here/ph_field.py"; verdict "FIELD CERT PASS"
if build -DDUMP "$here/cert_pk.c" -o "$tmp/dump" && "$tmp/dump" > "$tmp/dump.txt"; then check 1 python3 "$here/ph_oracle.py" "$tmp/dump.txt"; verdict "mismatches -> PASS"
else echo "FAIL building or running the field product dump"; fail=1; fi
for m in 8 64; do
  build -DPH_M=$m "$here/cert_pk.c" -o "$tmp/cpk" && check 1 "$tmp/cpk" || fail=1; verdict "CERT PK PASS"
  build -DPH_M=$m -DNEG "$here/cert_pk.c" -o "$tmp/cpk" || fail=1
  if "$tmp/cpk" > /dev/null; then echo "FAIL PK negative control (m=$m): a derivation with colliding exponents passed"; fail=1
  else echo "PASS PK negative control (m=$m): a derivation with colliding exponents is rejected"; fi
done
for m in 16 32 64; do build -DPH_M=$m "$here/cert_2l.c" -o "$tmp/c2l" && check 1 "$tmp/c2l" || fail=1; done
build -DNEG "$here/cert_2l.c" -o "$tmp/c2l" || fail=1
if "$tmp/c2l" > /dev/null; then echo "FAIL 2L negative control: a region function with colliding exponents passed"; fail=1
else echo "PASS 2L negative control: a region function with colliding exponents is rejected"; fi
check 1 python3 "$here/ph256_bounds.py"; verdict "score 255): PASS"
[ $fail = 0 ] && echo "ChainHash-256 certificates: ALL PASS" || { echo "ChainHash-256 certificates: FAIL"; exit 1; }
