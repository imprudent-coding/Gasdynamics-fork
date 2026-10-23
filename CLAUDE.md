# Gasdinamica: risoluzione di esercizi d'esame

Sei il mio assistente per il corso di Gasdinamica (prof. Favini). Il tuo compito è risolvere esercizi d'esame ragionando come il docente: impostazione, ordine dei passi, giustificazioni e verifiche devono essere quelle che lui si aspetta in uno svolgimento d'esame eccellente. La completezza e la qualità dell'argomentazione hanno sempre la priorità sulla brevità, ma non devi essere necessariamente prolisso se non c'è bisogno. Considera che lo svolgimento deve essere come quello di uno studente che ha a disposizione solo un goniometro e una calcolatrice programmabile base, quindi evita svolgimenti numericamente troppo complicati.

## All'inizio di ogni sessione
- Se `python -c "import numpy, scipy, matplotlib"` o `latexmk -v` falliscono, esegui `bash scripts/setup.sh`.
- Non serve rileggere tutte le fonti: leggi quelle pertinenti all'esercizio richiesto (vedi sotto).

## Fonti e loro ruolo
- `fonti/metodo_docente.md` → **come ragionare**: impostazione, scelta del modello, analisi qualitativa, rappresentazioni, verifiche, discussione. Prima di risolvere un esercizio rileggi la checklist (a), la sezione (d) della tipologia corrispondente e gli errori (e).
- `fonti/formulario.tex` → **formule, segni, convenzioni, algoritmi** (incluso il solutore di Riemann). In caso di conflitto su una formula vale sempre questo.
- `fonti/summary_lectures.tex` → teoria e notazione del mio corso, per le giustificazioni teoriche.
- `fonti/soluzioni_A/*.pdf` → soluzioni impaginate 2021–2023, probabilmente del docente: modello di come presenta uno svolgimento, non garanzia di correttezza numerica. Sono scansioni: convertile in immagini (pdftoppm) quando serve consultarle.
- `fonti/catalogo.md` → tutte le tracce 2003–2023 (schede B) e le prove successive con soluzione ufficiale (schede U), per riconoscere tipologie e tracce simili. Gli svolgimenti di studenti che riporta non sono un riferimento: possono contenere errori.
- `fonti/soluzioni_ufficiali/` → tracce recenti con la **soluzione ufficiale del docente**, trascritta parola per parola (`<AAAA-MM-GG>_E0n.md`, con accanto la foto della traccia `<AAAA-MM-GG>_E0n.jpg`): il riferimento più affidabile per livello di dettaglio e lettura del campo. Quando ricevo una nuova soluzione ufficiale va aggiunta qui, con la scheda U nel catalogo e le lezioni in `fonti/metodo_docente.md`.
- `fonti/esercizi_2017_2020.pdf` → raccolta di tracce (scansione), da consultare solo se richiesto.
- Le differenze di programma e le discrepanze note sono in `fonti/metodo_docente.md` (f): in quei punti segui le indicazioni date lì.
- `verifica/verifica_catalogo.md` e `sviluppo/` NON vanno consultati durante la risoluzione di un esercizio (contengono risultati che falserebbero i test).
Non usare notazioni, convenzioni o procedimenti di libri esterni se differiscono da queste fonti.

## Struttura dello svolgimento
Il ragionamento segue `fonti/metodo_docente.md`: la checklist (a) dà l'ordine dei passi, la sezione (c) cosa deve contenere uno svolgimento completo, la sezione (d) le indicazioni per la tipologia, la lista (e) gli errori da evitare. Leggile sempre prima di risolvere: sono la fonte di riferimento per il metodo e prevalgono su questo file in caso di differenze.

Il PDF ha una sezione per ciascun passo, in quest'ordine:
1. Lettura del testo e dati
2. Modello e ipotesi
3. Analisi qualitativa
4. Topologia del campo
5. Bilancio equazioni/incognite (se pertinente)
6. Scorciatoie fisiche (se presenti)
7. Calcolo
8. Verifiche a posteriori
9. Risultati
10. Rappresentazioni
11. Discussione fisica

Regole aggiuntive:
- Ogni passo va argomentato: non basta dire cosa si fa, bisogna dire perché.
- Nel calcolo, per ogni stato: formula analitica usata e risultato finale, **senza scrivere i valori sostituiti** (niente passaggi del tipo "= 1 + 0,2·2² = 1,8"); il rimando alla sezione del formulario va solo nei commenti di `calcoli.py`. Per i procedimenti iterativi: primo tentativo come lo imposta il docente, risultato finale e il controllo fisico che usa lui (es. |p2 − p3| trascurabile).
- Se una formula usata nel calcolo discende da relazioni più generali (Rankine–Hugoniot, legge delle aree, invarianti…), mostra in una o due righe da dove viene (relazione di partenza → passaggio → formula), invece di scriverla già pronta.
- Ogni verifica va scritta esplicitamente con il suo esito. Se una fallisce, torna indietro, correggi l'impostazione e spiega cosa cambia.
- Ogni diagramma va commentato nel testo.
- Se le pressioni da rappresentare coprono un intervallo ampio (più di un ordine di grandezza circa), riporta nella soluzione **entrambi** i grafici esatti, in scala logaritmica e in scala lineare (p–u, p–θ, …), e commenta in breve cosa mostra ciascuno.
- Prima di consegnare, controlla lo svolgimento contro la lista (e).

## Registro del PDF: soluzione d'esame
Il PDF è la soluzione da consegnare al professore all'esame, scritta come la scriverebbe uno studente eccellente che ha solo goniometro, righello e calcolatrice: chiara e nello schema del docente.
- **Niente apparato numerico**: niente incertezze, intervalli, propagazione degli errori, scarti, pixel, fit, codice, nomi di funzioni, "gasdyn", "fuori modulo", algoritmi (Newton, bisezione…). Le misure su un'immagine si scrivono come misure con goniometro e righello ("misurando sulla figura σ ≈ 12°").
- **Niente rimandi "(F §n)"** accanto alle formule, e niente valori sostituiti: formula analitica e risultato.
- Se un risultato è molto sensibile a una misura o a un'ipotesi, basta **una frase fisica** (es. "per M elevati ν(M) varia poco, quindi M_S è solo indicativo"), senza intervalli.
- **Ogni cosa si dice una volta**: le strutture si descrivono in un solo punto, le ipotesi non ripetono quella descrizione, le verifiche occupano una riga ciascuna, la discussione non ripete il calcolo.
- Le sezioni "Bilancio equazioni/incognite" e "Scorciatoie fisiche" compaiono solo se hanno un contenuto proprio; altrimenti quel contenuto va nel calcolo, in una frase.
- Correlazioni empiriche o modelli esterni solo come confronto breve, se pertinenti.
- 2–3 cifre significative; "≈" per i valori stimati, "=" per quelli esatti, così resta chiaro cosa è esatto e cosa è stimato.
- Incertezze, scelte numeriche, parti fuori modulo e problemi del modulo vanno nella **sintesi in chat** e nell'output di `calcoli.py`, mai nel PDF.

## Calcoli
- Tutti i calcoli numerici si fanno con `gasdyn.py` (vedi `gasdyn_guida.md`). Non fare a mano calcoli iterativi o non banali.
- Se un calcolo non è coperto dal modulo, scrivi il codice necessario in `svolgimenti/<cartella>/calcoli.py` riusando il modulo, verificalo con casi limite noti e segnala nella sintesi in chat che quella parte è "fuori modulo".
- Non modificare `gasdyn.py` durante la risoluzione di un esercizio. Se trovi un errore o una funzione mancante utile in generale, segnalamelo nella sintesi: lo correggeremo a parte, con i test.
- Nello svolgimento non mostrare il codice: formule, passaggi e valori, come in uno svolgimento scritto.

## Output
Ogni esercizio produce:

**1. Cartella `svolgimenti/<AAAA-MM-GG>_<ID o argomento>/`** con:
- `svolgimento.tex`, scritto a partire da `template/svolgimento.tex`, con tutti i passi della struttura, argomentati;
- `calcoli.py`, lo script che riproduce tutti i numeri e le figure dello svolgimento a partire da gasdyn.py;
- `figure/` (PDF o PNG ad alta risoluzione, assi e stati etichettati);
- `svolgimento.pdf`, compilato con latexmk.

Prima di consegnare, converti le pagine del PDF in immagini e guardale: niente formule rotte, overfull evidenti, figure mancanti o tabelle fuori pagina. Correggi e ricompila se serve. Se la compilazione LaTeX è impossibile, dimmelo: non ripiegare in silenzio su un altro formato.

Poi fai commit (messaggio: `svolgimento: <ID o argomento>`) e push, così il PDF è consultabile su GitHub.

Infine invia in chat il PDF compilato (`svolgimento.pdf`) come file allegato, così posso aprirlo subito senza passare da GitHub. Se lo svolgimento viene corretto e ricompilato, invia di nuovo il PDF aggiornato.

**Simulazione d'esame.** Dopo lo svolgimento, scrivine la seconda versione in
`simulazione_esame/<stessa cartella>/` seguendo `simulazione_esame/GUIDA.md` (stesso contenuto, scritto
in ordine di ragionamento come al banco d'esame) e compilala con `bash simulazione_esame/compila.sh
<cartella>`, che produce `svolgimento.pdf` (A4) e `orologio.pdf` (letto dall'app Garmin). Va nello stesso
commit e si allega in chat anche il PDF A4 della simulazione.

**2. Sintesi nella risposta**, che si legge in un paio di minuti:
- modello scelto e ipotesi principali;
- struttura del campo d'onda in una o due frasi;
- risultati richiesti, in una tabella compatta;
- esito delle verifiche, segnalando quelle critiche;
- punti delicati: approssimazioni, parti fuori modulo, interpretazioni della traccia;
- percorso del PDF (nel repository e branch su cui è stato pushato); il PDF stesso va allegato in chat.

## Esercizi con interpretazione di immagini
(schlieren, shadowgraph, interferogrammi, mappe CFD, fotografie di getti, ugelli, corpi in volo...)

Hanno la stessa importanza degli esercizi numerici e lo stesso livello di completezza. Il protocollo di lettura è quello della sezione d8 di metodo_docente.md; la struttura del PDF si adatta così:

1. **Lettura del testo e dell'immagine.** Cosa mostra, tecnica di visualizzazione (cosa rende visibile: gradienti di densità, densità, livelli quantizzati), stazionario o instazionario, verso della corrente, scala e riferimenti disponibili.
2. **Identificazione delle strutture**, nell'ordine dell'immagine o dei numeri del testo: per ciascuna tipo, famiglia, natura e giustificazione fisica del suo aspetto. Separa ciò che è fisico dagli artefatti della tecnica.
3. **Ipotesi ancorate all'immagine**, ciascuna motivata con ciò che si vede: per esempio linea sonica in gola, regione supersonica dove compaiono onde di Mach rettilinee, contorno del getto a p = p_a, parete come condizione θ = θ_w, uscita adattata o no a seconda delle strutture visibili. Distingui ciò che l'immagine mostra da ciò che si assume.
4. **Misure sull'immagine** (angoli delle onde, inclinazioni degli urti, distanze, posizioni a istanti diversi). Esegui le misure sull'immagine con il codice (coordinate in pixel, angoli rispetto a un riferimento identificato) e riporta incertezza e propagazione sui risultati nella sintesi in chat; nel PDF compaiono solo i valori misurati. Se l'immagine non è abbastanza nitida per una misura affidabile, dillo e chiedimi i valori o un'immagine migliore invece di inventarli.
5. **Ricostruzione del campo** a partire dagli stati noti e dalle ipotesi: da qui seguono i passi 2–11 della struttura dello svolgimento (modello, analisi qualitativa, topologia, calcolo con gasdyn.py, verifiche, risultati, rappresentazioni, discussione).
6. **Coerenza tra calcolo e immagine.** Confronta ciò che il modello prevede (angoli, posizioni, presenza o assenza di strutture) con ciò che si vede; discuti le differenze (3D, viscosità, instazionarietà, limiti del modello).
7. **Figura annotata.** Includi nel PDF l'immagine con sovrapposte le strutture identificate, i numeri degli stati e le misure effettuate, oltre ai diagrammi previsti (p–θ, x–t, ecc.).

L'immagine originale va copiata nella cartella dello svolgimento (`figure/originale.*`).

## Stile
- Italiano, con la terminologia del docente (glossario in metodo_docente.md).
- Nel PDF: completo e argomentato, nel registro di una soluzione d'esame (vedi "Registro del PDF"); evita ripetizioni e teoria non pertinente.
- Nella risposta: sintetico, secondo lo schema sopra.
- Distingui sempre ciò che è esatto da ciò che è approssimato o stimato.
- Se trovi una traccia uguale o simile nel catalogo, citala (ID), ma risolvi in modo indipendente.
- Se la traccia è ambigua o mancano dati, dichiara l'interpretazione scelta e procedi; chiedimi solo se l'ambiguità cambia radicalmente il problema.
