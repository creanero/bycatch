#!/bin/bash
# Runs on the compute host from a copy of this repo (see CLAUDE.local.md for paths).
# All WFI chips of both POTS-1 nights: field-only extraction, then a two-night period search per chip.
export PYTHONPATH=~/exotic-gpu/exotic-private-git OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=~/miniforge3/envs/exotic_gpu_env/bin/python
OUT=~/bycatch-data/results/eso/wfi/pots_field
mkdir -p $OUT
ln -sf ~/bycatch-data/results/check/pots-1_54173_float.json $OUT/54173_chip3.json
ln -sf ~/bycatch-data/results/check/pots-1_54230_float.json $OUT/54230_chip3.json
for chip in 1 2 4 5 6 7 8; do
  for night in 54173 54230; do
    f=$OUT/${night}_chip${chip}.json
    [ -s $f ] && continue
    echo "== extract night $night chip $chip $(date)"
    $PY -W ignore -m pipeline.wfi extract ~/bycatch-data/eso/wfi/${night}_OTSF-1A --chip $chip --workers 16 \
        --aperture-coeff 1.2 --out $f || echo "FAILED night $night chip $chip"
  done
done
for chip in 1 2 3 4 5 6 7 8; do
  echo "== periodsearch chip $chip $(date)"
  mkdir -p $OUT/search_chip$chip
  $PY -W ignore -m pipeline.periodsearch $OUT/54173_chip$chip.json $OUT/54230_chip$chip.json \
      --out $OUT/search_chip$chip/candidates.json --plot-dir $OUT/search_chip$chip --workers 16 || echo "FAILED search chip $chip"
done
echo "== done $(date)"
