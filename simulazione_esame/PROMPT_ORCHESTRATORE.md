# Prompt per l'orchestratore: soluzione completa di una prova d'esame

Da incollare in una sessione "orchestratore" che ha accesso a questa repository e a
`imprudent-coding/garmin-latex`, allegando le foto della prova. L'orchestratore:
1. salva le foto in `tracce/`, ritaglia un'immagine per esercizio e scrive l'indice;
2. crea un thread per esercizio con le istruzioni complete;
3. segue i thread e alla fine consegna in chat i PDF e una tabella riassuntiva.

Ogni thread fa tutto il flusso di `CLAUDE.md`:
- svolgimento di riferimento;
- simulazione d'esame (`GUIDA.md`);
- compilazione A4 + orologio;
- merge su `main`;
- verifica della pubblicazione sull'orologio.

**Unica riga da cambiare: `PROVA`.**
- Può essere una data vera (`del 14/09/2026`) o una data fittizia (`del 01/02/1970`): in questo caso usane una diversa per ogni prova, così sull'orologio restano in gruppi separati.
- Può essere anche un nome (`foto 3`), se garmin-latex supporta le prove senza data.

Se l'orchestratore non riesce a salvare le foto come file, si ferma al punto 1a: in quel
caso si allegano le foto direttamente nei thread.

```
PROVA: del 01/02/1970        ← unica riga da cambiare (data vera, fittizia, o un nome)
Allego le foto della prova (una o più; una foto può contenere più esercizi, e un
esercizio può continuare su due foto).

Organizza la soluzione completa della prova "Prova PROVA" con un thread per esercizio.
Tu non risolvi nulla: prepari, distribuisci e controlli.

1. Tracce nel repository (imprudent-coding/Gasdynamics-fork). <P> = data di PROVA in formato
   AAAA-MM-GG (se PROVA è un nome: prova-<nome>).
   a. Salva le foto come file in tracce/<P>/foto1.<ext>, foto2.<ext>, … nell'ordine in cui
      le ho allegate. Se non riesci a salvarle come file, fermati e dimmelo.
   b. Guarda le foto e individua gli esercizi, con la numerazione della prova. Per ognuno
      salva tracce/<P>/E0n.png: un ritaglio leggibile che contenga tutto il testo e tutte
      le figure di quell'esercizio e niente degli altri (con un po' di margine; se
      l'esercizio è su due foto, metti i pezzi uno sotto l'altro nello stesso file).
   c. Scrivi tracce/<P>/indice.md: per ogni esercizio, il numero, la foto o le foto di
      origine e la posizione ("foto1, metà superiore"), e una riga con l'argomento.
   d. Se qualcosa è ambiguo (numerazione poco chiara, esercizio tagliato, parti
      illeggibili), fermati e chiedimi prima di andare avanti.
   e. Commit "tracce: Prova PROVA", PR verso main e merge. Mostrami l'indice.

2. Un thread per esercizio. Dopo il merge, per ogni esercizio n crea una sessione sul
   repository imprudent-coding/Gasdynamics-fork (branch main), nello stesso ambiente di
   questa sessione, titolo "Prova PROVA – E0n", con questo messaggio iniziale (sostituisci
   PROVA, n e <P>):

   "Risolvi l'esercizio n della Prova PROVA. La traccia è tracce/<P>/E0n.png su main; le foto
   originali sono nella stessa cartella e tracce/<P>/indice.md dice da dove viene il
   ritaglio: se il ritaglio è tagliato o poco nitido, usa la foto originale.
   Segui CLAUDE.md integralmente. Ricava una volta:
   - CARTELLA = <P>_E0n_<argomento breve>;
   - TITOLO = \esercizio{Prova PROVA -- Esercizio n}{<argomento>}.
   0. Prepara l'ambiente come dice CLAUDE.md (setup.sh se serve) e leggi CLAUDE.md,
      simulazione_esame/GUIDA.md e, in fonti/metodo_docente.md, la checklist (a), la (c),
      la (d) della tipologia e la (e). Se la traccia non è leggibile, fermati e dillo.
   1. Svolgimento di riferimento in svolgimenti/CARTELLA/ (11 passi, gasdyn.py,
      calcoli.py, figure, la traccia copiata come figure/originale.*, PDF controllato
      pagina per pagina e contro la lista (e); cita tracce simili del catalogo).
   2. Simulazione d'esame in simulazione_esame/CARTELLA/ secondo GUIDA.md, con TITOLO
      esattamente così: leggibile dall'inizio alla fine senza perdere nulla (confronta
      sezione per sezione con lo svolgimento). Compila con compila.sh, guarda l'A4 e
      scorri orologio.pdf cercando righe fuori colonna.
   3. Commit "svolgimento: Prova PROVA E0n <argomento>", push, PR verso main e merge. Se
      main è avanzato per gli altri esercizi, portalo nel tuo branch prima.
   4. Verifica la pubblicazione: workflow esami-orologio della fork riuscito; ultimo
      workflow notes su main di imprudent-coding/garmin-latex riuscito (se annullato da
      uno più recente, segui quello); l'ultima Release notes-* contiene il gruppo della
      Prova PROVA con la voce En.
   5. Chiudi con la sintesi prevista da CLAUDE.md, l'esito della verifica e i percorsi dei
      due PDF A4 (svolgimento e simulazione)."

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
   Se un PDF manca su main (thread fallito o merge non fatto), dimmelo nella tabella
   invece di saltarlo in silenzio.

Non fare altre modifiche ai repository oltre a tracce/.
```

Sull'orologio, a pubblicazione finita: menu → *Mostra* → *Simulazioni d'esame*, poi
*Sincronizza ora*.
