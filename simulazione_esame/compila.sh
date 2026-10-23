#!/usr/bin/env bash
# Compila una simulazione d'esame nei due formati (vedi GUIDA.md):
#   svolgimento.pdf  A4
#   orologio.pdf     colonna stretta, per l'app dell'orologio
# Uso: bash simulazione_esame/compila.sh <cartella> [<cartella> ...]
#      (senza argomenti: tutte le cartelle)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
dirs=("$@")
if [ ${#dirs[@]} -eq 0 ]; then
  mapfile -t dirs < <(find . -mindepth 2 -maxdepth 2 -name svolgimento.tex -printf '%h\n' | sort)
fi
for d in "${dirs[@]}"; do
  d="${d%/}"
  echo "== $d"
  (
    cd "$here/$d"
    latexmk -pdf -interaction=nonstopmode -halt-on-error -silent svolgimento.tex
    latexmk -pdf -interaction=nonstopmode -halt-on-error -silent -jobname=orologio \
      -usepretex='\def\formatoorologio{}' svolgimento.tex
    latexmk -c -silent svolgimento.tex
    latexmk -c -silent -jobname=orologio svolgimento.tex
  )
done
