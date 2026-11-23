#!/usr/bin/env bash
# Run the built FlexAID binary on the tiny complex in this directory.
#
# The repository does not ship a docking example or a test suite. These
# files are a synthetic Ala-Gly chain plus a three-heavy-atom ligand,
# converted with ProcessLigand (atom index 90000, residue 9999). The run
# is one GA generation and must exit 0 and write a pose PDB.
#
# COMPLF VCT is the mode that allocates the contact buffers. Scoring
# always calls Vcontacts, so this is the path a real docking uses.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$script_dir/../.." && pwd)"
binary="${1:-$root/WRK/FlexAID}"

if [[ ! -x "$binary" ]]; then
  echo "FlexAID binary is not executable: $binary" >&2
  exit 1
fi

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/state" "$work/tmp"

cat > "$work/config.txt" << EOF
PDBNAM ${script_dir}/target.inp.pdb
INPLIG ${script_dir}/ligand.inp
METOPT GA
OPTIMZ 9999 - -1
OPTIMZ 9999 - 0
RNGOPT LOCCEN 2.500 1.500 0.500 6.000
IMATRX ${root}/WRK/MC_st0r5.2_6.dat
DEPSPA ${root}/WRK
MAXRES 1
COMPLF VCT
STATEP ${work}/state
TEMPOP ${work}/tmp
EOF

"$binary" "$work/config.txt" "$script_dir/ga.txt" "$work/smoke"

test -s "$work/smoke_INI.pdb"
test -s "$work/smoke_0.pdb"
grep -q HETATM "$work/smoke_0.pdb"

echo "smoke test passed: $(wc -c < "$work/smoke_0.pdb") bytes written to smoke_0.pdb"
