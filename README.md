# Gasdynamics_Favini

Ambiente per risolvere con Claude Code gli esercizi d'esame di Gasdinamica (prof. Favini):
fonti del corso, metodo del docente, modulo di calcolo validato `gasdyn.py` e template LaTeX
per gli svolgimenti. Le istruzioni per Claude sono in [`CLAUDE.md`](CLAUDE.md).

## Struttura

| Percorso | Contenuto |
| --- | --- |
| `CLAUDE.md` | Istruzioni di risoluzione: fonti, struttura dello svolgimento, output |
| `gasdyn.py` | Modulo di calcolo (numpy, scipy, matplotlib) con la notazione del corso |
| `gasdyn_guida.md` | Guida d'uso del modulo: convenzioni, funzioni, esempi, limiti |
| `requirements.txt` | Dipendenze Python |
| `fonti/` | Formulario, summary delle lezioni, metodo del docente, catalogo delle tracce 2003–2023, raccolta 2017–2020, `soluzioni_A/` (17 soluzioni 2021–2023) |
| `template/svolgimento.tex` | Template LaTeX degli svolgimenti |
| `simulazione_esame/` | Seconda versione di ogni svolgimento, scritta come al banco d'esame (A4 + formato orologio per l'app Garmin): vedi `GUIDA.md` |
| `svolgimenti/` | Una sottocartella per esercizio: `svolgimento.tex/.pdf`, `calcoli.py`, `figure/` |
| `scripts/setup.sh` | Installazione idempotente dell'ambiente (Ubuntu): Python, LaTeX, poppler |
| `tests/` | Test del modulo (livelli 1–3) |
| `verifica/` | Validazione del modulo sugli svolgimenti del corso (non consultare durante la risoluzione) |
| `sviluppo/` | Materiale di sviluppo: note, prompt, trascrizioni delle lezioni e appunti |

## Uso

```bash
bash scripts/setup.sh          # una volta per ambiente
python -m pytest -q tests      # 352 test
```

Poi, in una sessione di Claude Code: incollare la traccia (o indicarne l'ID del catalogo, o allegare
l'immagine). Lo svolgimento viene creato in `svolgimenti/<AAAA-MM-GG>_<ID>/`, compilato e pubblicato
con un commit.

`fonti/esercizi_2003_2023.pdf` (113 MB) non è nel repository: il suo contenuto è riassunto in
`fonti/catalogo.md`.
