# Catalogo preliminare per l'estensione multi-paziente

Data della ricognizione: 2026-10-23

## Scopo

Questo documento raccoglie coorti pubbliche candidate per una validazione
multi-paziente di InSituGame. La priorità è mantenere separati:

- paziente biologico;
- campione/blocco;
- sezione;
- regione;
- spot o cellula;
- piattaforma e batch.

Gli spot o le cellule aumentano l'informazione intra-campione, ma non
sostituiscono le repliche biologiche indipendenti. Nessuna coorte qui sotto è
stata ancora scaricata o analizzata localmente.

Questo non implica che lo studio originale sia invalido perché usa pochi
campioni. La numerosità adeguata dipende dalla domanda:

- per descrivere eterogeneità, composizione e relazioni spaziali **dentro una
  sezione**, migliaia di cellule o spot possono essere informativi;
- per stimare una proprietà che deve generalizzare a pazienti diversi, il
  paziente è il livello biologico rilevante e serve una valutazione separata;
- coorti esterne non sono automaticamente una validazione migliore: differenze
  di tessuto, istotipo, trattamento, fissazione, piattaforma, risoluzione e
  pipeline possono introdurre un effetto di dominio.

L'estensione viene quindi trattata come test di **trasferibilità e robustezza**,
non come correzione automatica dello studio originale.

## Candidati prioritari

| Dataset | Indicazione e numerosità verificata | Piattaforma | Utilità potenziale | Accesso | Stato |
|---|---|---|---|---|---|
| [E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530) | NSCLC: 8 pazienti con tumore; anche tessuto adiacente non coinvolto e 2 donatori sani. Il record descrive 20 sezioni tumorali, 16 adiacenti e 4 da donatori sani. | 10x Visium | Confronto tumore/tessuto adiacente/controlli e test di stabilità tra pazienti. | Record BioStudies/ArrayExpress pubblico; i raw e i processati sono collegati alla struttura EBI/ENA. | **Metadati verificati**; sample sheet ancora da normalizzare localmente. |
| [GSE300445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300445&targ=self&form=text&view=quick) | Melanoma cutaneo primario: 4 tumori/pazienti Visium. Lo studio include inoltre bulk RNA-seq da 57 campioni pretreatment e snRNA-seq da 2 metastasi. | 10x Visium, immunofluorescenza ad alta dimensionalità, bulk RNA-seq, snRNA-seq | La combinazione spaziale e ortogonale è adatta a confrontare interazioni tumore-immunità e annotazioni proteiche. | GEO pubblico con archivio supplementare indicato nel record. | **Candidato forte**, ma quattro pazienti sono ancora pochi per inferenze cliniche. |
| [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Atlante di tumori mammari; il catalogo GEO conferma la serie single-cell/spatial. La componente Visium è riportata come sei campioni tumorali in fonti collegate allo studio. | 10x Visium e single-cell/CITE-seq | Buon candidato per breast cancer e per definire firme cellulari/CAF prima dell'analisi spot-level. | GEO pubblico; i file e i metadati effettivamente disponibili devono essere controllati prima del download. | **Candidato utile**, da verificare campione per campione. |
| [GSE246518](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246518&targ=self&form=text&view=quick) | Tessuto mammario: 6 portatrici BRCA2 e 5 non portatrici, come riportato nel disegno della serie. | 10x Visium | Può fornire una coorte multi-campione per testare robustezza e confondimento genetico. | Processati e metadati pubblici; i file raw sono dichiarati mancanti per ragioni di privacy. | **Candidato condizionale**: non adatto se il protocollo richiede FASTQ raw. |
| [E-MTAB-14560](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14560) | Breast, lung e DLBCL; il record pubblico espone 20 campioni e 32 assay, ma non un conteggio univoco dei pazienti nella pagina principale. | Visium CytAssist, GeoMx DSP, Chromium Flex | Utile per validazione cross-platform e per separare effetti biologici da effetti di piattaforma. | BioStudies pubblico; dati di sequenziamento collegati a [ERP165490](https://www.ebi.ac.uk/ena/browser/view/ERP165490). | **Candidato metodologico**, solo dopo deduplicazione dei pazienti dal sample sheet. |

## Fonte di espansione

Il [portale HTAN](https://data.humantumoratlas.org/explore) è una fonte
prioritaria per ampliare la coorte. Il portale distingue HTAN Phase 1 e Phase
2 e usa modelli dati differenti. La documentazione ufficiale
([HTAN Manual](https://docs.humantumoratlas.org/)) indica che l'accesso e i
livelli di controllo dipendono dal progetto e dal tipo di file.

HTAN è promettente perché può collegare caso/paziente, biospecimen, assay
spaziali, single-cell, patologia e imaging. Tuttavia, in questa fase non è
ancora verificato:

- quali progetti contengano simultaneamente coordinate, counts raw,
  annotazioni di regione e marker ortogonali;
- quanti pazienti indipendenti siano disponibili per ciascun istotipo;
- se i file necessari siano open access o richiedano autorizzazioni;
- se i dati siano armonizzabili senza mescolare piattaforme incompatibili.

## Falso positivo escluso

La ricerca iniziale ha associato il PMID del lavoro sul melanoma spaziale a
[GSE113244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113244).
La verifica diretta tramite NCBI E-utilities mostra però che GSE113244 è una
serie di RNA non codificante su cellule di glioma U87, non la coorte di
melanoma spaziale. Non deve quindi essere usata come accession per questa
estensione.

Questo controllo è un esempio del motivo per cui ogni accession deve essere
verificata direttamente nel record ufficiale, non solo tramite risultati di
ricerca o riassunti generati automaticamente.

## Cosa è stato verificato

1. I record GEO/EBI indicati esistono e sono associati alle serie o agli studi
   descritti.
2. E-MTAB-13530 è una coorte NSCLC Visium multi-paziente con tessuto tumorale,
   adiacente e controlli come descritto dal record.
3. GSE300445 è una coorte melanoma con componente Visium e dati ortogonali
   riportati nel catalogo della ricerca.
4. GSE176078 è la serie GEO dell'atlante single-cell/spatial del breast
   cancer.
5. GSE246518 dichiara una componente raw non disponibile per motivi di privacy.
6. E-MTAB-14560 espone campioni e assay multi-piattaforma, ma non consente di
   dedurre in sicurezza il numero di pazienti dalla sola pagina principale.
7. HTAN separa Phase 1 e Phase 2 e documenta l'esistenza di requisiti di
   accesso dipendenti dai dati.
8. L'audit locale dei soli SDRF ha confermato 160 righe di assay e 40 sorgenti
   in E-MTAB-13530, con 10 identificativi individuali espliciti: 8 pazienti e
   2 donatori. Le righe assay sono quindi ripetizioni tecniche/sezioni e non
   160 pazienti.
9. L'SDRF di E-MTAB-14560 contiene 128 righe e 20 sorgenti, ma nessuna colonna
   esplicita `patient`, `individual`, `subject` o `donor`. I nomi delle sorgenti
   suggeriscono gruppi come B1, B2, L1 e L2, ma non sono sufficienti per
   dedurre un conteggio affidabile di pazienti.
10. Per GSE300445 il series matrix ufficiale conferma quattro campioni Visium,
    la designazione come tumori primari e la disponibilità dell'archivio raw;
    non è stata scaricata la taratura raw né alcuna immagine.

## Cosa è ancora da verificare localmente

L'audit preliminare è stato eseguito in
`/home/luca/Scrivania/ESTENSIONE_PILOT_METADATA/` e ha scaricato soltanto:

- `study.json` per E-MTAB-13530 e E-MTAB-14560;
- i due file SDRF;
- il series matrix di GSE300445 (circa 3 KB).

Il riepilogo strutturato è
`/home/luca/Scrivania/ESTENSIONE_PILOT_METADATA/metadata_audit.tsv`.
Questi file sono artefatti di lavoro separati e non fanno parte del repository
principale.

Per un pilot paired minimo sono stati scaricati da E-MTAB-13530 soltanto i
file spaziali e le matrici filtrate di `P10_B1` e `P10_T1`, entrambi dello
stesso paziente:

- `P10_B1`: tessuto adiacente, 1.755 barcode filtrati, 36.601 feature;
- `P10_T1`: tumore, 2.200 barcode filtrati, 36.601 feature;
- entrambi hanno 4.992 righe nella griglia spaziale e sovrapposizione completa
  dei barcode filtrati con le coordinate;
- mediane UMI/feature rilevate: 989/654 per `P10_B1` e 4.434/2.192 per
  `P10_T1`.

Questa differenza di profondità è un controllo QC, non una prova biologica.
Il pilot paired serve a verificare importazione, coordinate e confronto
intra-paziente; non consente di stimare generalizzazione clinica con un solo
paziente.

### Output esplorativo dell'assay

Applicando l'implementazione corrente InSituGame ai due campioni, senza
modificare il codice:

| Campione | Regime stimato | R² in-sample |
|---|---|---:|
| `P10_B1` adiacente | Resistant Dominance | 0.0151 |
| `P10_T1` tumore | Stag Hunt (Bistable), `x_sep = 0.206` | 0.0201 |

Questi numeri sono **solo un controllo esplorativo**, non un risultato
biologico. In particolare, l'R² molto basso indica che il modello corrente
spiega una frazione ridotta del fitness proxy in questi due campioni. Il cambio
di regime non può essere attribuito al tumore rispetto all'adiacente perché
sono presenti contemporaneamente differenze di profondità, composizione,
microambiente e campionamento. Servono replicati e analisi di sensibilità
prima di interpretarlo.

I risultati completi sono in
`/home/luca/Scrivania/ESTENSIONE_PILOT_METADATA/pilot_assay_results.json`.

Prima di qualsiasi analisi statistica bisogna scaricare solo i metadati e
verificare:

- identificativo paziente non ambiguo e assenza di duplicati;
- identificativo di sezione e di campione;
- coordinate spaziali;
- matrice counts e formato;
- annotazioni tumorali, stromali e immunitarie;
- qualità e soglie di filtraggio;
- batch, centro, piattaforma e preparazione;
- licenza e condizioni di riuso;
- disponibilità di marker T_R e CAF, oppure di firme cellulari validate;
- corrispondenza tra dati spaziali e dati ortogonali.

Il numero di spot dichiarato nei dataset Visium non deve essere usato come
numero di repliche. La validazione primaria dovrà usare split e bootstrap a
livello paziente; la sezione potrà essere un livello gerarchico intermedio.

## Disegno pilota raccomandato

### Fase 1: audit dei metadati

1. Partire da E-MTAB-13530 e GSE300445.
2. Scaricare soltanto sample sheet, metadata e liste file.
3. Costruire una tabella:
   `patient_id -> sample_id -> section_id -> region_id -> spot/cell_id`.
4. Escludere qualsiasi campione per cui il paziente non sia ricostruibile.
5. Prima di unire dataset diversi, stimare separatamente le quantità nel
   dataset originale e in quello esterno; non fare pooling automatico.

### Fase 2: validazione statistica

La prima analisi non deve sostituire il risultato originale con una coorte
eterogenea. Deve produrre tre risultati distinti:

1. **replicazione interna**: stabilità del metodo dentro la coorte originale;
2. **replicazione intra-dominio**: pazienti aggiuntivi con piattaforma e
   definizione biologica comparabili;
3. **trasferimento cross-domain**: comportamento su piattaforme o tessuti
   diversi, dichiarato esplicitamente come analisi di generalizzazione.

Confrontare, con split interamente a livello paziente:

- M1;
- M2c;
- M2d;
- M3;
- un modello gerarchico con effetto paziente;
- un modello con covariate di piattaforma/batch quando necessario.

Riportare performance per paziente, intervalli bootstrap a livello paziente,
eterogeneità tra pazienti e confronto con null model spaziali. Non usare un
unico pooled fit sugli spot come evidenza di generalizzazione.

### Fase 3: controllo dei confondenti

- ripetere l'analisi con firme T_R/CAF alternative;
- separare composizione cellulare da interazione spaziale;
- controllare correlazioni tra T_R e CAF;
- includere densità cellulare, regione patologica e batch;
- eseguire permutazioni che preservino la struttura spaziale;
- per l'optimal transport, costruire null che conservino la dipendenza
  comune da `x_S`.

### Fase 4: validazione ortogonale

Usare immunofluorescenza, CITE-seq, snRNA-seq o annotazioni patologiche solo
come validazione concordante. Anche una concordanza alta resta osservazionale:
non dimostra causalità né, da sola, che i coefficienti siano payoff biologici.

## Decisione operativa

Il primo pilot dovrebbe usare:

1. **E-MTAB-13530** per il polmone, perché offre più pazienti e controlli
   tissutali;
2. **GSE300445** per il melanoma, perché offre dati spaziali e ortogonali;
3. **GSE176078** come supporto breast/spatial e riferimento single-cell;
4. **HTAN** solo dopo aver stabilito un filtro di inclusione riproducibile.

Questa catalogazione non modifica il codice, il manoscritto, il repository
GitHub o il form di submission. Non è ancora una dimostrazione che una delle
coorti sia pronta per l'analisi: è un elenco verificato di candidati e dei
controlli che precedono il download e il modelling.

## Test basati su dati pubblicati

È legittimo usare accessioni pubbliche della letteratura come fixture di test,
ma con tre livelli distinti:

1. **test tecnico**: apertura di HDF5/MTX, allineamento barcode-coordinate,
   conteggio feature e conservazione delle masse;
2. **test di riproducibilità**: riesecuzione della pipeline su uno stesso
   campione con output deterministico;
3. **test biologico esterno**: replicazione o trasferibilità, da dichiarare
   come analisi separata e non come prova causale.

Le fixture candidate sono:

- E-MTAB-13530, già usato per il pilot paired;
- GSE300445, melanoma Visium con dati ortogonali;
- GSE176078, atlante breast cancer single-cell/spatial.

Il repository principale supera attualmente i 19 test `unittest`; il runner
`pytest` non è installato nell'ambiente globale e non è stato forzato con
installazioni invasive. I test sono stati eseguiti con:

`python3 -m unittest discover -s tests -v`

Tutti i 19 test sono passati. Questo conferma la coerenza del codice esistente
e dei benchmark matematici, non la validità scientifica delle stime sul nuovo
pilot.
