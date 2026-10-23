# Simulazioni d'esame

Questa cartella contiene, per ogni esercizio, una **seconda versione** dello svolgimento di
`svolgimenti/`: stessi numeri, stesse figure, stessa correttezza, ma scritta come la si scriverebbe
al banco d'esame, da leggere dall'inizio alla fine senza tornare indietro. Serve a vedere *come
formulare* la soluzione. Ogni esercizio si compila in due formati: A4 e "orologio", che l'app
Garmin (repository `garmin-latex`) mostra nel gruppo *Simulazioni d'esame*.

## Quando si scrive

Dopo lo svolgimento in `svolgimenti/<AAAA-MM-GG>_<ID>/` (che resta quello di riferimento,
nella struttura in 11 passi di `CLAUDE.md`), e con lo stesso nome di cartella:

```
simulazione_esame/<AAAA-MM-GG>_<ID>/svolgimento.tex   sorgente
simulazione_esame/<AAAA-MM-GG>_<ID>/svolgimento.pdf   A4
simulazione_esame/<AAAA-MM-GG>_<ID>/orologio.pdf      colonna stretta, per l'orologio
```

Compilazione (entrambi i formati): `bash simulazione_esame/compila.sh <AAAA-MM-GG>_<ID>`.
Per risolvere una prova intera con un thread per esercizio (dalle foto della traccia fino
alla pubblicazione sull'orologio) c'è il prompt pronto in `PROMPT_ORCHESTRATORE.md`.
I numeri vengono dal `calcoli.py` dello svolgimento di riferimento e le figure dalla sua
cartella `figure/` (con `\graphicspath`): non si ricalcola e non si duplica nulla.

## Come si scrive

Il contenuto è quello dello svolgimento di riferimento; cambiano ordine e tono.

- **Nessuna perdita di dettaglio.** La simulazione non è un riassunto: deve contenere tutto ciò
  che c'è nello svolgimento di riferimento, cioè ogni stato calcolato, formula, derivazione in una
  o due righe, ipotesi con la sua giustificazione, misura, verifica con il suo esito, figura e
  diagramma con il commento, e ogni punto della discussione fisica. Più leggibile non vuol dire
  più corta: si risparmia solo sulle ripetizioni che la struttura in 11 passi impone (la stessa
  struttura descritta in più sezioni), mai su un passaggio o un approfondimento.
  **Prima di consegnare**, rileggi lo svolgimento di riferimento sezione per sezione e controlla
  che ogni elemento compaia nella simulazione; se manca, aggiungilo dove serve nel ragionamento.

- **Ordine del ragionamento, non delle categorie.** Si segue il campo nell'ordine in cui lo si
  affronta: cosa si vede o cosa dice il testo → modello → poi, una regione o un evento alla volta
  (nell'ordine fisico: dalla gola a valle, dall'istante iniziale in avanti), *perché* succede, la
  formula, il risultato, il controllo. Non ci sono sezioni separate per topologia, verifiche,
  risultati e rappresentazioni.
- **Le verifiche stanno dove servono**: subito dopo il valore che verificano ("controllo subito
  che l'urto sia sul ramo debole…"), con l'esito.
- **I diagrammi stanno dove si usano**, commentati lì (x–t quando si descrivono le onde, p–u o
  p–θ quando si risolve l'interazione, il grafico richiesto alla fine).
- **Tono da compito scritto**: prima persona ("misuro", "assumo", "controllo"), frasi chiare,
  ogni passaggio motivato con un "perché" o un "quindi". Il lettore è il docente.
- **Chiusura**: la risposta alla domanda (grafico o tabella dei risultati) e la
  discussione fisica completa dello svolgimento di riferimento, senza ripetere il calcolo.
- Valgono tutte le regole del **registro del PDF** di `CLAUDE.md`: niente apparato numerico, niente
  valori sostituiti (formula e risultato), 2–3 cifre significative, "≈" per le stime, la
  derivazione in una o due righe delle formule che vengono da relazioni più generali, ogni cosa
  detta una volta. Come sempre, controllare contro la lista (e) di `fonti/metodo_docente.md`.
- La traccia va in testa, breve. Le sezioni (`\section*`) hanno titoli che dicono cosa si fa
  ("A parete: il Mach alla separazione"), non il nome del passo.

## Regole del sorgente (per il formato orologio)

Il preambolo comune è `preambolo.tex` (`\input{../preambolo}`). Nel formato orologio il testo sta in
una colonna di 4,6 cm: quello che la supera viene rimpicciolito sull'orologio, con lo zoom.

- Titolo con `\esercizio{Prova del GG/MM/AAAA -- Esercizio n}{argomento}`: va anche nel titolo del
  PDF, che l'orologio mostra nell'elenco.
- Figure con `\figura[frazione della larghezza A4]{file}{didascalia\label{...}}` (sempre `[H]`):
  nel formato orologio diventano una pagina con lo zoom.
- Più relazioni nella stessa formula: `\begin{gather*} A \ee B \ee C \end{gather*}`. In A4 stanno su
  una riga separate da virgole, sull'orologio una sotto l'altra. Evitare `\qquad` fra relazioni.
- Formule lunghe: spezzarle (risultato su una riga a parte con `\ee`); tabelle con poche colonne.
- Niente `\newpage`, colonne affiancate o `minipage` per affiancare testo e figure.

Prima di consegnare: guardare le pagine di `svolgimento.pdf` (come per gli svolgimenti) e scorrere
`orologio.pdf` cercando righe che escono dalla colonna senza motivo.
