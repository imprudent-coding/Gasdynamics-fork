# Prompt per l'orchestratore: soluzione completa di una prova d'esame

Due prompt da dare, uno dopo l'altro, alla stessa sessione "orchestratore", che deve avere
accesso a questa repository e a `imprudent-coding/garmin-latex`.

1. **Preparazione.** Senza foto. L'orchestratore crea un thread per esercizio. Ogni thread
   prepara l'ambiente, carica il contesto delle repository e riceve subito tutte le istruzioni
   di lavoro; poi aspetta la traccia.
2. **Tracce.** Con le foto allegate. L'orchestratore salva le foto in `tracce/`, ritaglia
   un'immagine per esercizio, la pubblica su `main` e dice a ogni thread di partire. Poi segue
   i thread e alla fine consegna in chat i PDF e una tabella riassuntiva.

Le istruzioni complete, merge compreso, stanno nel messaggio iniziale di ogni thread (prompt 1).
I messaggi che l'orchestratore manda dopo a un thread arrivano come messaggi di un'altra
sessione, non tuoi: per questo il prompt 2 dice solo "la traccia è pronta, procedi".

**In entrambi i prompt si cambiano solo le prime due righe, con gli stessi valori.**
- `PROVA`:
  - una data vera (`del 14/09/2026`);
  - una data fittizia (`del 01/02/1970`), diversa per ogni prova, così sull'orologio restano
    in gruppi separati;
  - oppure un nome (`foto 3`), se garmin-latex supporta le prove senza data.
- `ESERCIZI`: il numero di esercizi della prova.

## Prompt 1: preparazione dei thread

```
PROVA: del 01/02/1970        ← data vera, fittizia, o un nome
ESERCIZI: 3                  ← numero di esercizi della prova

Prepara la soluzione della prova d'esame di Gasdinamica "Prova PROVA": un thread di lavoro
per ciascuno degli ESERCIZI esercizi. Tu non risolvi nulla. Le foto delle tracce te le do
con il prossimo messaggio: per ora non serve nessuna traccia.

<P> = data di PROVA in formato AAAA-MM-GG (se PROVA è un nome: prova-<nome>).

Per ogni esercizio n = 1…ESERCIZI crea una sessione sul repository
imprudent-coding/Gasdynamics-fork (branch main), nello stesso ambiente di questa sessione,
titolo "Prova PROVA – E0n", con questo messaggio iniziale (sostituisci PROVA, n e <P>):

  "Questo thread risolve l'esercizio n della Prova PROVA. La traccia arriverà su main come
  tracce/<P>/E0n.png: quando un messaggio ti dice che è pronta, fai git pull di main e
  procedi. Fino ad allora fai solo la preparazione, senza modificare file.

  PREPARAZIONE (subito):
  - esegui le verifiche di inizio sessione di CLAUDE.md e, se servono, bash scripts/setup.sh;
    controlla python (numpy, scipy, matplotlib), latexmk e pdftoppm;
  - prova bash simulazione_esame/compila.sh su una delle cartelle già presenti in
    simulazione_esame/, poi ripristina i file (git checkout -- .);
  - leggi CLAUDE.md, simulazione_esame/GUIDA.md e, in fonti/metodo_docente.md, la checklist
    (a) e le sezioni (c) ed (e);
  - controlla di poter leggere imprudent-coding/garmin-latex.
  Poi rispondi «Pronto per l'esercizio n» con eventuali problemi, e aspetta.

  LAVORO (quando la traccia è pronta): segui CLAUDE.md integralmente. La traccia è
  tracce/<P>/E0n.png; le foto originali sono nella stessa cartella e tracce/<P>/indice.md
  dice da dove viene il ritaglio: se è tagliato o poco nitido usa la foto originale. Se la
  traccia non è leggibile, fermati e dillo. Ricava una volta:
  - CARTELLA = <P>_E0n_<argomento breve>;
  - TITOLO = \esercizio{Prova PROVA -- Esercizio n}{<argomento>}.
  1. Rileggi in fonti/metodo_docente.md la sezione (d) della tipologia. Svolgimento di
     riferimento in svolgimenti/CARTELLA/ (11 passi, gasdyn.py, calcoli.py, figure, la
     traccia copiata come figure/originale.*, PDF controllato pagina per pagina e contro la
     lista (e); cita tracce simili del catalogo).
  2. Simulazione d'esame in simulazione_esame/CARTELLA/ secondo GUIDA.md, con TITOLO
     esattamente così: leggibile dall'inizio alla fine senza perdere nulla (confronta
     sezione per sezione con lo svolgimento). Compila con compila.sh, guarda l'A4 e scorri
     orologio.pdf cercando righe fuori colonna.
  3. Commit "svolgimento: Prova PROVA E0n <argomento>", push, PR verso main e merge. Se
     main è avanzato per gli altri esercizi, portalo nel tuo branch prima.
  4. Verifica la pubblicazione: workflow esami-orologio della fork riuscito; ultimo
     workflow notes su main di imprudent-coding/garmin-latex riuscito; l'ultima Release
     notes-* contiene il gruppo della Prova PROVA con la voce En.
  5. Chiudi con la sintesi prevista da CLAUDE.md, l'esito della verifica e i percorsi dei
     due PDF A4 (svolgimento e simulazione)."

Quando tutti hanno risposto (controlla tu lo stato delle sessioni: non possono scriverti),
dammi una tabella: esercizio, link della sessione, stato ("pronto" o il problema).
Ricorda gli identificativi delle sessioni: ti servono per il prossimo messaggio.
Non creare altre sessioni e non fare commit.
```

## Prompt 2: tracce e avvio (allegando le foto)

```
PROVA: del 01/02/1970        ← gli stessi valori del prompt 1
ESERCIZI: 3

Allego le foto della Prova PROVA (una o più; una foto può contenere più esercizi, e un
esercizio può continuare su due foto). Sono ESERCIZI esercizi, già assegnati ai thread
che hai preparato.

1. Tracce nel repository (imprudent-coding/Gasdynamics-fork). <P> = data di PROVA in formato
   AAAA-MM-GG (se PROVA è un nome: prova-<nome>).
   a. Salva le foto come file in tracce/<P>/foto1.<ext>, foto2.<ext>, … nell'ordine in cui
      le ho allegate. Se non riesci a salvarle come file, fermati e dimmelo.
   b. Individua gli esercizi con la numerazione della prova. Per ognuno salva
      tracce/<P>/E0n.png: un ritaglio leggibile che contenga tutto il testo e tutte le
      figure di quell'esercizio e niente degli altri (con un po' di margine; se è su due
      foto, metti i pezzi uno sotto l'altro nello stesso file).
   c. Scrivi tracce/<P>/indice.md: per ogni esercizio, il numero, la foto o le foto di
      origine e la posizione ("foto1, metà superiore"), e una riga con l'argomento.
   d. Se gli esercizi trovati non sono ESERCIZI, o qualcosa è ambiguo (numerazione poco
      chiara, esercizio tagliato, parti illeggibili), fermati e chiedimi.
   e. Commit "tracce: Prova PROVA", PR verso main e merge. Mostrami l'indice.

2. Avvio. Dopo il merge, manda a ogni thread "Prova PROVA – E0n" il messaggio: "La traccia
   tracce/<P>/E0n.png è su main: fai git pull e procedi con il LAVORO del tuo messaggio
   iniziale."

3. Controllo e consegna. Segui i thread fino alla fine (stato della sessione ed eventi):
   non possono scriverti, quindi controlla tu periodicamente. Se uno si ferma con una
   domanda o un errore, riportamelo subito. Quando hanno finito tutti:
   a. aggiorna il tuo checkout di main (git pull) e allegami in chat, come file, i PDF di
      ogni esercizio: svolgimenti/<cartella>/svolgimento.pdf e
      simulazione_esame/<cartella>/svolgimento.pdf, in ordine di esercizio, con una
      didascalia "E0n – svolgimento" / "E0n – simulazione d'esame";
   b. dammi una tabella: esercizio, argomento, link della sessione, PR unita,
      pubblicazione sull'orologio (sì/no), e per ciascuno i risultati principali in una
      riga.
   Se un PDF manca su main (thread fallito o merge non fatto), dimmelo nella tabella invece
   di saltarlo in silenzio.

Non fare altre modifiche ai repository oltre a tracce/.
```

Sull'orologio, a pubblicazione finita: menu → *Mostra* → *Simulazioni d'esame*, poi
*Sincronizza ora*.

Se l'orchestratore non riesce a salvare le foto come file (prompt 2, punto 1a), alleghi le foto
direttamente in ciascun thread, scrivendo: "Ecco la traccia: procedi con il LAVORO del tuo
messaggio iniziale, usando la foto allegata al posto di tracce/<P>/E0n.png".
