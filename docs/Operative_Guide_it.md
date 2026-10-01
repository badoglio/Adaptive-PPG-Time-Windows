# Guida operativa alla dashboard

[English](Operative_Guide_en.md) · [Italiano](Operative_Guide_it.md) · [日本語](Operative_Guide_ja.md)

Questa guida spiega ogni controllo, grafico e tabella della dashboard web di **adaptive_ppg 0.2.0**. Per ogni
elemento dice che cosa significa, come sceglierne il valore per un dato studio e come leggere i risultati in
termini fisiologici. Per installare e avviare la dashboard, vedi la
[guida all'installazione](INSTALL_it.md). Il metodo è riassunto nel [README](../README.md).

* La dashboard è uno strumento di ricerca. Non è un dispositivo medico, e i suoi risultati non sono validati
  per la diagnosi.
* La guida descrive ciò che fa il codice della versione 0.2.0. Dove il README differisce dal codice, la guida
  segue il codice e lo segnala ([sezione 20](#20-problemi-noti)).
* Le etichette della GUI sono scritte in **grassetto**, esattamente come appaiono sullo schermo (in inglese).
  Le chiavi di configurazione (per il file JSON e la riga di comando) sono scritte come `codice`.
* Un **battito** è un singolo impulso, dal suo piede (onset) al piede dell'impulso successivo. Un **battito
  valido** è un battito che ha superato tutti i controlli di qualità. **u.a.** significa unità arbitrarie,
  cioè le unità del segnale in ingresso. **DS** significa deviazione standard.
* I numeri usano il punto decimale, come nella GUI e nei file esportati.

---

## Indice

1. [Come funziona l'analisi](#1-come-funziona-lanalisi)
2. [Flusso di lavoro consigliato](#2-flusso-di-lavoro-consigliato)
3. [Barra laterale: 1. Data](#3-barra-laterale-1-data)
4. [Barra laterale: 2. Method](#4-barra-laterale-2-method)
5. [Barra laterale: 3. Windows per category](#5-barra-laterale-3-windows-per-category)
6. [Barra laterale: 4. Signal processing](#6-barra-laterale-4-signal-processing)
7. [Barra laterale: 5. Fixed-window baselines](#7-barra-laterale-5-fixed-window-baselines)
8. [Intestazione, avvisi e KPI](#8-intestazione-avvisi-e-kpi)
9. [Scheda: Signal & beats](#9-scheda-signal--beats)
10. [Scheda: Variability & windows](#10-scheda-variability--windows)
11. [Scheda: Features](#11-scheda-features)
12. [Scheda: Templates](#12-scheda-templates)
13. [Scheda: Benchmark](#13-scheda-benchmark)
14. [Scheda: Export](#14-scheda-export)
15. [Guida alle feature](#15-guida-alle-feature)
16. [Impostazioni per caso d'uso](#16-impostazioni-per-caso-duso)
17. [Lista dei controlli](#17-lista-dei-controlli)
18. [Cosa riportare](#18-cosa-riportare)
19. [Valori di riferimento dei dati di esempio](#19-valori-di-riferimento-dei-dati-di-esempio)
20. [Problemi noti](#20-problemi-noti)
21. [Bibliografia](#21-bibliografia)

---

## 1. Come funziona l'analisi

La dashboard esegue la stessa pipeline del comando `adaptive-ppg run`.

1. **Caricamento.** Il PPG (e, se presente, un ECG) viene letto dal file. I valori mancanti (NaN, Inf) sono
   riempiti per interpolazione lineare. Le registrazioni più brevi di 10 s vengono rifiutate.
2. **Saturazione (clipping).** Alla frequenza nativa, i tratti di almeno 40 ms (e di almeno 2 campioni) che
   restano entro l'1 % dell'escursione del segnale dal minimo o dal massimo globale sono marcati come saturati.
3. **Ricampionamento** a una frequenza di elaborazione comune (125 Hz per default), con un filtro polifase.
4. **Filtraggio** in due rami, entrambi a fase zero (Butterworth di ordine 4 applicato in avanti e
   all'indietro al segnale meno la sua mediana):
   * un ramo di *rilevamento*, 0.5–8 Hz, usato solo per trovare i picchi sistolici;
   * un ramo di *morfologia*, 0.5–30 Hz per default, usato per tutto il resto.
5. **Polarità.** Se la derivata del segnale di rilevamento ha asimmetria (skewness) negativa, il segnale viene
   invertito (vedi **Signal polarity** nella [sezione 6](#6-barra-laterale-4-signal-processing)).
6. **Rilevamento dei picchi** (Elgendi et al. 2013 per default) e **segmentazione in battiti**. Ogni battito
   va dal suo onset all'onset del battito successivo, quindi l'ultimo picco rilevato, che non ha un onset
   successivo, viene scartato.
7. **Qualità dei battiti.** Ogni battito è confrontato con un riferimento costruito dai 30 battiti precedenti
   ed è marcato come valido o non valido.
8. **Variabilità morfologica.** Ogni battito viene normalizzato: si sottrae la retta dal suo onset alla sua
   fine, si scala il picco a 1 e si ricampiona il battito su 128 punti, con il picco sistolico al 30 % del
   battito. Una metrica di variabilità è calcolata sugli ultimi battiti validi normalizzati e trasformata in
   un indice **H** tra 0 (morfologia stabile) e 1 (morfologia che cambia).
9. **Finestre adattive.** Per ogni categoria di feature, H stabilisce la lunghezza della finestra (in battiti
   validi) e la sovrapposizione tra finestre consecutive: finestre lunghe quando la morfologia è stabile,
   finestre brevi quando cambia.
10. **Template d'insieme e feature.** I battiti di ogni finestra sono allineati e mediati in un template, e le
    feature sono misurate sul template. Le stesse feature sono misurate anche su ogni singolo battito e su
    finestre di durata fissa, per confronto.

### Equazioni

Metrica di variabilità di default, `drift_z`, calcolata sui battiti validi normalizzati (`N` = **N_past**):

```
drift_k   = RMS sui 128 campioni di [ media(battiti k−N+1 … k) − media(battiti k−2N+1 … k−N) ]
drift_z_k = drift_k / sqrt( (2/N) · σ²_k )
σ²_k      = varianza entro i blocchi dei 2N battiti, combinata sui due blocchi (2N − 2 gradi di libertà)
            e mediata sui 128 campioni
```

Se la morfologia è stazionaria e le fluttuazioni da battito a battito sono indipendenti, il valore atteso di
`drift²` è `2σ²/N`, quindi `drift_z ≈ 1` qualunque sia il livello di rumore. Un cambiamento a gradino di
ampiezza RMS Δ posto tra i due blocchi dà approssimativamente `drift_z ≈ sqrt(1 + N·Δ²/(2σ²))`. La metrica
misura quindi il cambiamento in unità della dispersione da battito a battito, e la sua sensibilità cresce con
√N.

Normalizzazione, lunghezza della finestra e sovrapposizione (lo stesso H guida tutte le categorie; ogni
categoria `c` ha i propri `T_min`, `T_max` e `T_crit`):

```
H      = clip( (v − lo) / (hi − lo), 0, 1 )           v = variabilità; lo, hi fissati da H normalization
s(H)   = 1 / (1 + exp( k · (H − θ) ))                  k = Sigmoid slope, θ = soglia
O      = O_max − (O_max − O_min) · s(H)                O_min = 0.25, O_max = 0.85
T_max  = max( T_max_sec / IBI_local , T_min )          in battiti; IBI_local = mediana delle durate degli ultimi 10 battiti validi
T*     = T_min + (T_max − T_min) · s(H)
L_k    = min( T*_k , L_(k−1) + Δ )                     crescita limitata a Δ battiti per battito (Δ = 1); la riduzione è immediata
T_k    = max( round(L_k), T_crit )                     lunghezza della finestra in battiti validi
step_k = max( 1, round( T_k · (1 − O_k) ) )            battiti tra due ancore di finestra consecutive
```

Valori del fattore di forma con i default (k = 10, θ = 0.3):

| H | 0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| s(H) | 0.953 | 0.881 | 0.731 | 0.500 | 0.269 | 0.119 | 0.001 |

Esempio svolto (`Sample1.CSV`, 71 bpm, IBI ≈ 0.84 s, categoria macro). `T_max` = 30 s / 0.84 s ≈ 35.5
battiti. Con H = 0, T ≈ 3 + 32.5 × 0.953 ≈ 34 battiti e O ≈ 0.28, quindi una nuova finestra inizia circa ogni
25 battiti. Con H = 1 la finestra si riduce subito a 3 battiti con O = 0.85, e una nuova finestra inizia a ogni
battito. Dopo il cambiamento, la finestra torna a crescere di al più un battito per battito (circa 31 battiti
per tornare da 3 a 34).

### Cosa H vede e cosa non vede

Poiché ogni battito è normalizzato prima di calcolare la variabilità:

* H **non** risponde a un puro cambiamento di ampiezza dell'impulso (tutti i campioni del battito scalati
  dello stesso fattore).
* H **non** risponde a un allungamento uniforme del battito nel tempo, e il ricampionamento a tratti (picco
  sempre al 30 %) elimina anche un cambiamento della proporzione tra salita e discesa. I soli cambiamenti di
  crest time e di duty cycle, quindi, non riducono le finestre.
* H **risponde** ai cambiamenti di forma all'interno della salita e all'interno della discesa: per esempio, un
  cambiamento dell'altezza o della posizione dell'onda tardo-sistolica/diastolica e dell'incisura dicrota
  rispetto alla discesa, o della curvatura della salita.
* Un cambiamento della frequenza cardiaca può comunque cambiare H indirettamente, perché la parte diastolica
  dell'impulso non si scala con la durata del ciclo.

Conseguenza pratica: quando la risposta di interesse è soprattutto un cambiamento di ampiezza o di tempi (per
esempio la riduzione dell'ampiezza dell'impulso con la vasocostrizione simpatica), H può restare basso e le
finestre restano lunghe. In questi protocolli, imposta un **T_max (s)** più breve per le categorie di interesse
([sezione 16](#16-impostazioni-per-caso-duso)).

---

## 2. Flusso di lavoro consigliato

1. **Carica** la registrazione ([sezione 3](#3-barra-laterale-1-data)). Leggi la riga informativa sotto il
   titolo: canale, frequenza di campionamento, durata, ECG, polarità.
2. Leggi gli **avvisi** (riquadri gialli). Se un avviso segnala una discrepanza della frequenza di
   campionamento o una colonna sbagliata, correggila prima di proseguire.
3. In **Signal & beats**, ingrandisci alcuni battiti. L'impulso filtrato deve puntare verso l'alto, con una
   salita ripida. I punti verdi devono stare sui picchi sistolici e i triangoli arancioni sui piedi. Gli
   intervalli rossi sono battiti scartati.
4. Controlla i quattro **KPI**: battiti validi, quota di battiti con H > θ, finestra mediana e SQI del
   template.
5. In **Variability & windows**, controlla che H salga dove ti aspetti un cambiamento (per esempio all'inizio
   di un compito) e resti basso a riposo.
6. In **Features**, guarda le feature di interesse con il loro CI95 e confrontale con le finestre fisse di
   riferimento.
7. In **Templates**, controlla che gli impulsi medi abbiano una forma plausibile.
8. In **Export**, scarica lo ZIP e il JSON di configurazione. Annota le finestre fisse e il preset: non sono
   salvati nel JSON.

Ogni modifica di un controllo riesegue l'analisi. I risultati sono tenuti in cache, quindi tornare a
un'impostazione precedente è rapido.

---

## 3. Barra laterale: 1. Data

**Source** sceglie da dove viene il segnale:

* **Synthetic demo (known ground truth)**: una registrazione simulata di cui si conoscono le feature vere.
  Serve alla scheda **Benchmark**.
* **Sample data**: le registrazioni della cartella `sample_data`. È la scelta di default quando la cartella
  esiste.
* **Upload a file**: una tua registrazione.

### Synthetic demo

* **Scenario**:
  * `stress`: 120 s di riposo, poi 120 s di stress acuto raggiunto rapidamente (costante di tempo 4 s), poi
    120 s di recupero lento (costante di tempo 25 s);
  * `rest`: 300 s, stazionario.
* **Sampling rate (Hz)**: 25–1000 Hz, default 125. Usalo per vedere come la frequenza di acquisizione
  influisce sulle feature, in particolare su quelle della derivata (SDPTG).
* **Noise seed**: 0–10000, default 1. Cambia il rumore casuale; la fisiologia sottostante resta la stessa.

Ogni battito simulato è la somma di tre onde gaussiane: sistolica, tardo-sistolica e diastolica. Il generatore
cambia questi parametri tra riposo e stress (ampiezze relative all'onda sistolica, ritardi dopo l'onda
sistolica):

| Parametro | Riposo | Stress |
| :--- | :---: | :---: |
| Frequenza cardiaca (bpm) | 70 | 100 |
| Onda sistolica: tempo / DS (s) | 0.16 / 0.055 | 0.13 / 0.045 |
| Onda tardo-sistolica: ampiezza / ritardo (s) / DS (s) | 0.35 / 0.09 / 0.045 | 0.12 / 0.08 / 0.045 |
| Onda diastolica: ampiezza / ritardo (s) / DS (s) | 0.45 / 0.24 / 0.09 | 0.22 / 0.19 / 0.08 |
| Tempo di arrivo dell'impulso (s) | 0.25 | 0.19 |
| Ampiezza dell'impulso | 1 | 0.6 |

Rumore e fattori di confondimento: rumore bianco con SNR di 30 dB; deriva della linea di base (0.3 × ampiezza
dell'impulso, 0.05–0.15 Hz); respirazione a 0.25 Hz con modulazione d'ampiezza dell'8 % e 3 bpm di aritmia
sinusale respiratoria; variazione casuale (jitter) dell'1 % dell'intervallo tra battiti; artefatti da
movimento circa una volta al minuto (1–4 s, 2 × ampiezza dell'impulso). È incluso un ECG sintetico, quindi il
PAT è disponibile. Il livello DC è 2.0, quindi l'indice di perfusione sintetico (circa 50 %) non è
fisiologico.

### Sample data

**File** elenca i file in `sample_data`, ordinati per nome. `Sample1.CSV` è il default.

| File | Fonte | Contenuto |
| :--- | :--- | :--- |
| `Sample1.CSV`, `Sample2.CSV` | Pulse Transit Time PPG Dataset 1.1.0 (Mehrgardt et al. 2022), licenza ODbL | 500 Hz. PPG infrarosso alla falange distale, soggetto seduto a riposo. Conteggi grezzi del MAX30101, che sono invertiti (più sangue, meno luce). Separatore `;`, virgola decimale. Circa il 23 % dei campioni consecutivi sono valori ripetuti. |
| `bidmc_01_Signals.csv` | Dataset BIDMC (Pimentel et al. 2017), derivato da MIMIC-II (Goldberger et al. 2000), licenza ODC-By | 125 Hz, 8 min, paziente critico. Colonne `Time [s]`, `RESP`, `PLETH`, `V`, `AVR`, `II`. La dashboard usa `PLETH` come PPG e la derivazione `II` come ECG. |

### Upload a file

**PPG recording** accetta file CSV, TXT, TSV, MAT, EDF e BDF fino a 200 MB.

* **CSV, TXT, TSV.** Il separatore (tabulazione, `;`, `,`, `|`, altrimenti spazi) e la virgola decimale sono
  rilevati automaticamente. Se il primo campo della prima riga è un numero, il file è letto senza intestazione
  e le colonne sono chiamate `col_0`, `col_1`, … Una colonna del tempo è riconosciuta dal nome: `time`, `t`,
  `sec`, `secs`, `seconds`, `tempo` o `timestamp` (maiuscole o minuscole), seguito da un confine di parola o
  da `_`, come in `Time [s]` o `time_ms`. I millisecondi sono riconosciuti da `[ms]`, `(ms)`, `_ms`, `msec` o
  `millisec` nel nome.
* **MAT.** Un nome sceglie la variabile che contiene il PPG. Un numero sceglie il canale della prima variabile
  numerica (per array 2-D). Una variabile chiamata `time` o `t` (o con un altro nome di tempo), altrimenti una
  chiamata `fs`, `sampling_rate`, `srate` o `sample_rate`, fornisce la frequenza di campionamento. Il campo
  dell'ECG è ignorato per i file MAT.
* **EDF, BDF.** La frequenza di campionamento è letta dall'intestazione. Il canale ECG deve avere la stessa
  frequenza del PPG. La lettura richiede il pacchetto `pyedflib` (opzione di installazione `edf`) o, in sua
  assenza, `mne`.

**PPG column / channel** — il nome della colonna o un indice a partire da 0. Se il campo è vuoto, si usa la
prima colonna diversa dal tempo (o il canale 0). I nomi sono confrontati prima esattamente, poi senza
distinguere maiuscole e minuscole; un numero è inteso come indice su tutte le colonne, compresa quella del
tempo. Per `bidmc_01_Signals.csv` il default è `PLETH`.

**ECG column / channel (optional, enables PAT)** — se indicato, i picchi R sono rilevati su questo canale e
viene calcolato il tempo di arrivo dell'impulso ([sezione 15](#15-guida-alle-feature)). Il rilevatore dei
picchi R è semplice: passa-banda 5–20 Hz, pendenza al quadrato smussata su 80 ms, un'unica soglia al 30 % del
suo 99° percentile per tutta la registrazione, picchi distanti almeno 0.3 s. È adeguato per derivazioni di
monitoraggio pulite, come la II di BIDMC. Per ECG ambulatoriali o rumorosi, controlla con cura i valori di PAT.
Per `bidmc_01_Signals.csv` il default è `II`.

**Infer the sampling rate from the file** — se selezionato, la frequenza viene dalla colonna del tempo o, per
EDF/BDF, dall'intestazione (sempre usata per EDF/BDF). La frequenza è stimata ai minimi quadrati su tutta la
colonna del tempo, quindi anche tempi arrotondati a 2–3 decimali danno la frequenza giusta. Se non è
selezionato, scrivi la frequenza in **Sampling rate (Hz)** (1–10000, default 125).

> **Suggerimento.** Una frequenza di campionamento sbagliata scala tutte le feature temporali e la frequenza
> cardiaca. Se la frequenza che scrivi differisce di più dell'1 % da quella implicata dalla colonna del tempo,
> un avviso indica il fattore di scala. Fidati della colonna del tempo, a meno che tu non sappia che è
> sbagliata.

---

## 4. Barra laterale: 2. Method

### Preset

**Preset** carica un insieme di valori di partenza. I controlli sottostanti li mostrano e possono modificarli.

| Preset | Differenze rispetto a `default` | Uso previsto |
| :--- | :--- | :--- |
| `default` | — | Uso generale |
| `rest` | θ = 0.5; limite di crescita Δ = 2 battiti per battito; **T_max (s)** 60 / 60 / 90 (macro / time_volume / derivatives) | Registrazioni a riposo: finestre lunghe che reagiscono solo a grandi cambiamenti |
| `acute_stress` | θ = 0.2; Δ = 0.5; **T_max (s)** 20 / 20 / 40 | Protocolli con transitori autonomici rapidi |
| `wearable_64hz` | **Morphology low-pass (Hz)** = 15; smussamento della derivata 110 ms; **T_crit** di derivatives = 15 | Dispositivi indossabili da polso campionati a circa 64 Hz |
| `legacy` | **Variability metric** = successive; **H normalization** = calibration; scala fissa 0.02–0.15; θ = 0.5 | Il metodo originale, per confronto |

> **Importante.** Il limite di crescita Δ e la finestra di smussamento della derivata non hanno un controllo
> nella GUI, quindi si applica sempre il valore del preset. I controlli di **3. Windows per category** **non**
> seguono un cambio di preset dopo la prima esecuzione: scrivi a mano i valori del preset per **T_max (s)** e
> **T_crit** ([sezione 20](#20-problemi-noti)). Gli altri controlli seguono il preset.

### Variability metric

| Metrica | Definizione | Comportamento |
| :--- | :--- | :--- |
| `drift_z` (default) | `drift` diviso per il suo valore atteso in condizioni stazionarie | ≈ 1 su un segnale stazionario, qualunque sia il livello di rumore; pensata per la scala fissa 1–3 |
| `drift` | Differenza RMS tra la media degli ultimi N battiti e la media degli N battiti precedenti | Misura il cambiamento persistente, in unità dell'impulso normalizzato (frazioni del picco); dipende dal livello di rumore |
| `successive` | Media retrospettiva su N battiti della differenza RMS tra battiti consecutivi | Misura soprattutto il rumore da battito a battito; usata dal preset `legacy` |
| `dispersion` | Deviazione RMS media degli ultimi N battiti dalla loro media | Dispersione entro il blocco, non cambiamento nel tempo |

`drift` e `drift_z` non sono definite finché non sono disponibili 2N battiti validi (nel frattempo H è tenuto
a 0). `successive` e `dispersion` sono definite dal secondo battito valido.

> **Suggerimento.** `drift`, `successive` e `dispersion` sono espresse come frazioni del picco (tipicamente
> qualche centesimo). Con **H normalization** = fixed e la scala di default 1–3, H resta a 0 e le finestre non
> si riducono mai. Usale con la normalizzazione calibration, rolling o global, oppure con una scala fissa dello
> stesso ordine (il preset `legacy` usa 0.02 e 0.15).

### H normalization

* **fixed** (default): H = 0 al valore **H = 0 at** e H = 1 al valore **H = 1 at**. La scala è la stessa per
  ogni registrazione, quindi H si può confrontare tra registrazioni e tra soggetti. Non c'è altro periodo di
  avvio (warm-up) oltre a quello della metrica stessa (i primi 2N − 1 battiti validi per `drift_z`).
* **calibration**: lo e hi sono il 5° e il 95° percentile della metrica durante i primi 60 s. Si assume che i
  primi 60 s siano una linea di base rappresentativa. Durante questi 60 s H è tenuto a 0 (warm-up). Se sono
  disponibili meno di 3 valori, si usa la scala fissa.
* **rolling**: lo e hi sono il 5° e il 95° percentile degli ultimi 120 battiti validi. La scala segue i
  cambiamenti lenti. Durante i primi 60 s H è tenuto a 0.
* **global**: percentili dell'intera registrazione. Usa dati futuri, quindi è solo per l'analisi offline.

> **Suggerimento.** Le normalizzazioni relative (calibration, rolling, global) stendono su 0–1 qualunque
> variabilità abbia la registrazione. Per costruzione, circa il 5 % dei battiti di riferimento raggiunge H = 1
> e una quota considerevole supera θ, anche quando non cambia nulla. Preferisci **fixed** con `drift_z`, a meno
> che tu non abbia un motivo per riscalare ogni registrazione.

**H = 0 at** e **H = 1 at** (mostrati solo con la normalizzazione fixed; tre decimali; **H = 1 at** è
mantenuto sopra **H = 0 at**). Con `drift_z`, il default 1 significa "nessun cambiamento oltre quello atteso
dalla dispersione da battito a battito", e 3 significa un cambiamento circa tre volte più grande. Sulle
registrazioni di esempio a riposo il `drift_z` mediano è 0.83–1.18
([sezione 19](#19-valori-di-riferimento-dei-dati-di-esempio)). Se una registrazione a riposo mostra molte
riduzioni delle finestre, guarda i suoi valori di `variability_used` nella **Beat quality table** e valuta di
alzare **H = 0 at** leggermente sopra la loro mediana.

### θ e pendenza

**θ (H at which windows start shrinking)** — 0.05–0.95, default 0.3. Con H = θ la lunghezza della finestra è a
metà tra `T_max` e `T_min`. Il valore corrispondente della metrica è `lo + θ·(hi − lo)`: con `drift_z` e la
scala 1–3, θ = 0.3 corrisponde a `drift_z` = 1.6. Un θ più basso fa ridurre le finestre per cambiamenti più
piccoli (più reattive, ma feature più rumorose). Un θ più alto le fa ridurre solo per grandi cambiamenti.

**Sigmoid slope k** — 2–30, default 10. Stabilisce quanto è netta la transizione tra finestre lunghe e brevi.
Con k = 10 e θ = 0.3, s(0) = 0.953, quindi le finestre sono quasi a `T_max` quando non cambia nulla. Con una
pendenza bassa la transizione è graduale, ma le finestre non raggiungono mai `T_max` (k = 2 dà s(0) = 0.646).
Con k ≥ 20 il comportamento è quasi binario (lunghe o brevi).

### N_past

**N_past (beats per drift block)** — 4–30, default 10. Il numero di battiti in ciascuno dei due blocchi
confrontati da `drift` e `drift_z`. È anche la lunghezza della media di `successive`, di `dispersion` e del
regressore IBI. La scelta è un compromesso:

* **N più grande**: più battiti in ogni media, quindi si rilevano cambiamenti persistenti più piccoli (la
  risposta di `drift_z` cresce come √N) e le modulazioni periodiche si mediano meglio.
* **N più piccolo**: risposta più rapida e warm-up più breve. La metrica richiede 2N battiti validi: con
  N = 10 il warm-up è di 19 battiti validi, circa 16 s a 70 bpm. Per un cambiamento a gradino, il `drift`
  retrospettivo raggiunge il massimo circa N battiti dopo il gradino e resta elevato per circa 2N battiti.

Indicazione fisiologica: ogni blocco dovrebbe coprire almeno circa **due cicli respiratori**, in modo che la
modulazione respiratoria dell'impulso si medi in entrambi i blocchi. In battiti,
`N ≥ 2 × frequenza cardiaca / frequenza respiratoria`.

* A 15 atti/min e 70 bpm un ciclo respiratorio dura circa 4.7 battiti, quindi N = 10 copre circa due cicli.
* Con respirazione lenta o guidata (6 atti/min a 60 bpm, 10 battiti per ciclo), usa N ≥ 20.
* Le onde di Mayer (circa 0.1 Hz; Julien 2006) durano circa 12 battiti a 70 bpm e con l'N di default si mediano
  solo in parte.

### Window anchor

**Window anchor** — `trailing` (default) o `centered`.

* **trailing**: ogni finestra termina al suo battito di ancoraggio. Ogni valore al tempo t dipende solo dai
  battiti fino a t (causale, come in tempo reale). Una media retrospettiva ritarda rispetto a un trend di
  circa metà della sua lunghezza, quindi i valori delle feature ritardano rispetto alla fisiologia di circa T/2
  (circa 15 s per una finestra di 30 s a riposo).
* **centered**: ogni finestra è centrata sul suo battito di ancoraggio. Il valore della feature non è ritardato
  rispetto all'ancora, il che è preferibile per l'analisi offline legata a eventi. H e l'IBI locale sono letti
  N_past//2 battiti più avanti, ma H reagisce comunque circa N/2 battiti dopo un cambiamento (invece di circa
  N battiti). Le finestre centrate usano dati futuri, quindi sono solo per l'analisi offline.

### Regress out IBI-driven variability

**Regress out IBI-driven variability** — se selezionato, viene rimossa la parte della metrica che segue le
fluttuazioni degli intervalli tra battiti. Il regressore è la media retrospettiva su N battiti di
`|ΔIBI| / IBI`. Una stima causale ai minimi quadrati, su una finestra che si allarga, inizia dopo 30 coppie, e
la metrica diventa `v − slope · (x − mean(x))`. La pendenza è salvata nei metadati (`ibi_regression_slope`).

Usalo quando le fluttuazioni degli intervalli tra battiti, per esempio una forte aritmia sinusale respiratoria
o un ritmo irregolare, fanno sì che H segua gli intervalli anziché la forma. Lascialo disattivato quando ti
aspetti un cambiamento reale di forma insieme a un cambiamento degli intervalli, perché parte del cambiamento
reale verrebbe rimossa.

---

## 5. Barra laterale: 3. Windows per category

Le feature sono raggruppate in tre categorie. Ogni categoria ha i propri limiti di finestra:

| Categoria | Feature | **T_min (beats)** | **T_max (s)** | **T_crit** |
| :--- | :--- | :---: | :---: | :---: |
| **macro** (17) | ampiezza dell'impulso, larghezze dell'impulso, aree, rapporti tra aree, ampiezza dell'inflessione, indice di riflessione, pendenze, livello DC, indice di perfusione | 3 | 30 | 3 |
| **time_volume** (10) | durata dell'impulso, frequenza cardiaca, crest time, tempo di discesa, duty cycle, tempo dell'incisura, ΔT_DVP, indice di rigidità, PAT al piede e al picco | 3 | 30 | 3 |
| **derivatives** (6) | SDPTG b/a, c/a, d/a, e/a, indice di invecchiamento, (b − e)/a | 5 | 60 | 10 |

* **T_min (beats)** (1–60): lunghezza della finestra, in battiti validi, quando la morfologia cambia (H vicino
  a 1). Stabilisce la risoluzione temporale durante i transitori: 3 battiti sono circa 2.5 s a 70 bpm.
* **T_max (s)** (5–300, passo 5): lunghezza della finestra, in secondi, quando la morfologia è stabile. È
  convertita in battiti con l'IBI locale (mediana degli ultimi 10 battiti validi), quindi a 100 bpm una
  finestra di 30 s contiene circa 50 battiti. Stabilisce lo smussamento a riposo e il ritardo delle finestre
  trailing (circa T_max/2).
* **T_crit** (1–60): il numero minimo di battiti validi in una finestra. Le finestre con meno battiti non sono
  prodotte, e la lunghezza della finestra non scende mai sotto T_crit.

Indicazioni:

* Le feature della derivata hanno finestre più lunghe perché la derivata seconda amplifica il rumore.
* Per i punti fiduciali fragili (incisura dicrota, rapporto delle aree al punto di inflessione, onde c e d
  della SDPTG), alza il **T_crit** della categoria.
* Scegli un **T_max (s)** più breve del cambiamento più breve che vuoi seguire senza che H se ne accorga. Questo
  conta per i cambiamenti di sola ampiezza o di soli tempi, che H non vede
  ([sezione 1](#cosa-h-vede-e-cosa-non-vede)).
* All'inizio della registrazione una finestra trailing contiene tutti i battiti validi disponibili fino a quel
  momento, quindi `n_valid` può essere minore di `T_beats`.
* Il limite di crescita Δ (`max_expand_per_beat`) non ha un controllo nella GUI: 1 per default, 2 con `rest`,
  0.5 con `acute_stress`. Per impostare un altro valore, usa un file di configurazione con la riga di comando.

---

## 6. Barra laterale: 4. Signal processing

**Resample to a common rate** (attivo per default) e **Processing rate (Hz)** (50–1000, passo 25,
default 125) — l'intera analisi avviene a questa frequenza. I template sono sempre costruiti su una griglia a
500 Hz (interpolazione B-spline cubica), quindi il filtro della derivata è lo stesso a ogni frequenza. Il
sovracampionamento non aggiunge banda: le frequenze di taglio dei filtri sono limitate al 45 % della minore tra
la frequenza di acquisizione e quella di elaborazione. Se l'opzione non è selezionata, si usa la frequenza
nativa (più lento per registrazioni a 500–1000 Hz).

> **Suggerimento.** Il README riporta, su `Sample1.CSV` e `Sample2.CSV`, una differenza mediana rispetto
> all'analisi nativa a 500 Hz, espressa in DS dei valori per battito, di 0.02 a 250 Hz, 0.10 a 125 Hz, 0.14 a
> 100 Hz e 0.3 a 64 Hz. 125 Hz è adeguato per la maggior parte delle feature. Per le feature della derivata,
> usa 250 Hz o più quando la frequenza di acquisizione lo consente. Frequenze sotto 100 Hz non sono consigliate
> per le feature della derivata.

**Signal polarity** — `auto` (default), `normal` o `inverted`. Sia nel PPG in trasmissione sia in quello in
riflessione, più sangue assorbe più luce, quindi la luce rilevata diminuisce durante la sistole. I conteggi
grezzi del fotorivelatore sono quindi capovolti rispetto al volume ematico (come in `Sample1.CSV` e
`Sample2.CSV`). I monitor clinici di solito mostrano la curva pletismografica già invertita, con il volume
verso l'alto. `auto` inverte il segnale quando la derivata del segnale di rilevamento ha asimmetria negativa:
in un impulso orientato correttamente la salita è più breve e più ripida della discesa, quindi la derivata ha
asimmetria positiva. Scegli `normal` o `inverted` quando conosci la convenzione, o quando `auto` sbaglia
(impulsi smorzati con salita e discesa quasi simmetriche, registrazioni molto rumorose).

**Morphology high-pass (Hz)** — 0.05–5, default 0.5. Rimuove la componente continua, la deriva della linea di
base e l'oscillazione respiratoria della linea di base. Poiché il filtro è applicato in avanti e all'indietro,
l'attenuazione alla frequenza di taglio è −6 dB.

| Frequenza (Hz) | 0.3 | 0.5 | 0.667 (40 bpm) | 0.75 | 1 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Attenuazione con passa-alto a 0.5 Hz | −36 dB | −6.0 dB | −0.76 dB | −0.29 dB | −0.02 dB |

* A 0.5 Hz, una frequenza del polso di 40 bpm perde 0.76 dB alla fondamentale, e una di 30 bpm perde 6 dB. Per
  la bradicardia (sotto circa 45 bpm), abbassa il passa-alto a 0.3 Hz (−0.13 dB a 0.5 Hz).
* Un passa-alto vicino alla frequenza del polso cambia la forma dell'impulso. Evita valori sopra circa 0.7 Hz.
* Valori più bassi lasciano più deriva della linea di base. La linea di base lineare da onset a fine, rimossa
  da ogni battito, ne compensa una parte.

**Morphology low-pass (Hz)** — 5–60, default 30 (−0.33 dB a 22.5 Hz, −6 dB a 30 Hz). La derivata seconda
amplifica le alte frequenze, quindi le onde della SDPTG dipendono da questo valore insieme allo smussamento
della derivata (70 ms) e alla frequenza di elaborazione. Usa 15 Hz per i dispositivi indossabili a bassa
frequenza di campionamento (preset `wearable_64hz`; −0.56 dB a 11.25 Hz). Valori più bassi smussano la salita e
cambiano i rapporti della SDPTG. Mantieni lo stesso valore all'interno di uno studio e riportalo. Se la
frequenza di taglio supera il 45 % della frequenza di campionamento, viene abbassata automaticamente e compare
un avviso.

**Bypass the morphology filter** — il ramo di morfologia diventa il segnale meno la sua mediana. Il ramo di
rilevamento resta filtrato. Usalo solo per ingressi già limitati in banda e privi di deriva della linea di
base, o per vedere l'effetto del filtro. Senza filtro, la deriva della linea di base entra nei template (in
parte rimossa dalla linea di base lineare di ogni battito) e le feature della derivata sono più rumorose.

**Onset method** — `tangent` (default) o `minimum`.

* **tangent**: l'intersezione della tangente nel punto di massima pendenza della salita con la retta
  orizzontale passante per il minimo che precede la salita (precisione inferiore al campione). Il piede a
  tangenti intersecanti è meno sensibile alla forma del piede e al rumore, ed è una scelta comune per i tempi
  di transito e di arrivo dell'impulso (Mukkamala et al. 2015).
* **minimum**: il minimo tra picchi consecutivi, raffinato con un'interpolazione parabolica.

Il piede a tangente di solito cade dopo il minimo. Con `tangent`, il crest time è più breve e `pat_foot` più
lungo che con `minimum`. Non mescolare i due metodi all'interno di uno studio.

**Min. correlation with the reference beat** — 0.5–0.99, default 0.86. Un battito supera `corr_ok` se la
correlazione di Pearson tra la sua forma normalizzata e il riferimento (mediana dei 30 battiti precedenti) è
almeno pari a questo valore. 0.86 riprende la soglia di template matching per il PPG di Orphanidou et al.
(2015), che la applicavano alla correlazione media su segmenti di 10 s; qui è applicata a ogni battito.
Li & Clifford (2012) discutono altri indici di qualità basati sul template matching per segnali pulsatili.

* Una soglia più alta scarta più battiti atipici (artefatti, battiti ectopici). Può anche scartare battiti
  durante un cambiamento di forma reale e rapido, il che sbilancia l'analisi verso la stazionarietà.
* Una soglia più bassa conserva più battiti ma lascia entrare artefatti nei template.

Gli altri controlli di qualità sono fissi: frequenza del polso 30–220 bpm; durata del battito entro ±30 % della
durata di riferimento; ampiezza tra 1/3 e 3 volte l'ampiezza di riferimento; nessuna saturazione; valori
finiti. Per i primi 5 battiti il riferimento è la mediana dei primi 30 battiti.

**Peak detector** — `elgendi` (default) o `biosppy`. `elgendi` è il rilevatore a due medie mobili di Elgendi
et al. (2013), con finestre di 111 ms e 667 ms e un offset di 0.02, applicato al ramo di rilevamento
0.5–8 Hz. La sua distanza minima tra picchi di 0.3 s lo limita a 200 bpm. `biosppy` esegue l'implementazione
BioSPPy dello stesso algoritmo sullo stesso ramo; richiede il pacchetto opzionale `biosppy` (e `peakutils`).
Con entrambi i rilevatori i picchi sono raffinati sul segnale di morfologia entro ±50 ms, con interpolazione
parabolica. La scheda **Benchmark** usa sempre `elgendi`.

**Subject height (m, 0 = unknown)** — abilita l'indice di rigidità, `SI = altezza / ΔT_DVP`
(Millasseau et al. 2002). L'altezza è usata come approssimazione della lunghezza del percorso arterioso.

---

## 7. Barra laterale: 5. Fixed-window baselines

**Fixed windows (s, comma-separated)** — default `10, 30`. Ogni valore crea una sequenza di finestre di durata
fissa, calcolate per tutte le 33 feature. Sono accettati sia le virgole sia i punti e virgola. Un campo vuoto
significa nessuna finestra fissa. Un testo che non si riesce a leggere mostra "Use numbers separated by
commas." e si usano 30 s.

Le finestre fisse usano lo stesso **Window anchor** delle finestre adattive:

* ancore ogni `W · (1 − overlap)` secondi, a partire dal primo picco più W (trailing) o più W/2 (centrata);
  l'ultimo picco è sempre un'ancora;
* una finestra contiene i battiti validi il cui picco cade in `(t_anchor − W, t_anchor]` (trailing) o entro
  `t_anchor ± W/2` (centrata);
* le finestre con meno di 3 battiti validi sono saltate.

**Overlap of fixed windows** — 0–0.9, default 0.5. Una sovrapposizione maggiore dà valori più frequenti ma lo
stesso smussamento.

Usa le finestre fisse come riferimento. Le finestre fisse brevi (10 s) seguono i cambiamenti rapidi ma sono
rumorose. Quelle lunghe (30–60 s) sono lisce ma ritardano e sfumano i transitori. Le finestre adattive cercano
di combinare i due comportamenti.

---

## 8. Intestazione, avvisi e KPI

### Intestazione e riga informativa

Il titolo "Adaptive-window PPG morphology" è seguito dalla versione del pacchetto. La riga informativa riporta
il nome del file, il canale PPG, la frequenza di acquisizione e, tra parentesi, la frequenza di elaborazione,
la durata, e le note ", with ECG (PAT enabled)" e ", inverted" quando si applicano.

### Avvisi ed errori

| Messaggio | Significato | Cosa fare |
| :--- | :--- | :--- |
| "Sampling rate X Hz is low: derivative features (SDPTG) will be unreliable." | Acquisizione sotto 50 Hz | Non interpretare le feature SDPTG |
| "Detection / Morphology high cutoff X Hz is too close to Nyquist …; using Y Hz." | Una frequenza di taglio del passa-basso superava il 45 % della frequenza di campionamento ed è stata abbassata | Riporta la frequenza di taglio effettivamente usata (`effective_highcut_hz` nei metadati) |
| "Signal inverted automatically (derivative skewness S < 0)." | La polarità `auto` ha invertito il segnale | Controlla in **Signal & beats** che l'impulso ora punti verso l'alto |
| "Only X% of beats passed quality checks." | Meno del 50 % di battiti validi | Controlla la colonna, la frequenza, la polarità e la qualità del segnale |
| "Sampling rate set to X Hz but the time column implies Y Hz. …" | La frequenza scritta e la colonna del tempo differiscono di più dell'1 % | Usa la frequenza del file, a meno che la colonna del tempo non sia sbagliata |
| "No beats detected." | Nessun picco trovato | Colonna sbagliata, segnale piatto o frequenza sbagliata |
| "Fewer than 20 valid beats: check the column, the sampling rate and the signal polarity." | Troppo pochi battiti da analizzare (errore) | Come sopra |
| "Could not analyse the recording: …" | La pipeline si è fermata, per esempio "Signal too short (x s); at least 10 s are required." | Leggi il messaggio |
| "Invalid configuration: …" | Impostazioni incoerenti, per esempio passa-alto ≥ passa-basso | Correggi l'impostazione indicata nel messaggio |

### KPI

* **Valid beats** — battiti validi / battiti rilevati, con la percentuale. Le registrazioni di esempio danno
  98–100 %. Sotto circa il 90 %, guarda gli intervalli scartati prima di interpretare le feature.
* **Beats with H > θ** — la quota di battiti, dopo il warm-up, in cui H supera θ, cioè in cui le finestre sono
  più brevi del punto intermedio tra `T_max` e `T_min`. I battiti non validi sono contati con l'H dell'ultimo
  battito valido. Con le impostazioni di default, le registrazioni di esempio a riposo danno 7–26 % e la
  registrazione sintetica di stress (seed 1) dà 14.8 %. Un valore alto a riposo suggerisce una reale non
  stazionarietà (postura, schema respiratorio, onde vasomotorie lente, contatto del sensore) oppure un θ o un
  **H = 0 at** troppo bassi.
* **Median macro window** — la lunghezza mediana della finestra macro su tutti i battiti (warm-up compreso), in
  battiti, e in secondi (lunghezza mediana × IBI locale mediano). Confrontala con **T_max (s)**.
* **Median template SQI** — la mediana, sulle finestre macro, della correlazione di Pearson media tra ogni
  battito e il template della finestra. Le registrazioni pulite danno valori vicini a 1 (0.985–0.993 sui dati
  di esempio). Ogni battito contribuisce al template con cui è confrontato, quindi l'SQI è distorto verso
  l'alto nelle finestre brevi.

---

## 9. Scheda: Signal & beats

**Show the raw signal** (attivo per default) aggiunge un secondo pannello con il segnale grezzo.

Il grafico mostra:

* **Filtered (morphology branch)**: il segnale usato per tutte le feature, dopo la correzione della polarità.
  Il titolo riporta la banda di morfologia, la banda di rilevamento e la frequenza di elaborazione
  ("PPG — unfiltered" quando il filtro è escluso).
* **Raw**: il segnale ricampionato nella sua polarità originale. **Non** è invertito, quindi con i conteggi
  grezzi punta verso il basso.
* Punti verdi: picchi sistolici dei battiti validi. Punti rossi: picchi dei battiti scartati. Triangoli
  arancioni: onset.
* Intervalli rossi ombreggiati: battiti scartati. Intervalli viola ("artefact"): artefatti simulati (solo dati
  sintetici).
* Le linee con più di 20 000 punti sono decimate per la visualizzazione (min–max). I marcatori sono esatti.

### Beat quality table

Una riga per ogni battito rilevato.

| Colonna | Significato |
| :--- | :--- |
| `onset_idx`, `peak_idx`, `end_idx` | Onset, picco sistolico e fine (= onset successivo), come indici di campione alla frequenza di elaborazione |
| `t_onset`, `t_peak` | Tempi dell'onset e del picco (s) |
| `duration_s` | Durata del battito, dall'onset all'onset successivo (s): l'intervallo del polso |
| `bpm` | 60 / `duration_s` |
| `amplitude` | Segnale di morfologia al picco meno quello all'onset (u.a.) |
| `template_corr` | Correlazione di Pearson del battito normalizzato con il battito di riferimento |
| `amplitude_ratio` | Ampiezza divisa per l'ampiezza di riferimento |
| `finite_ok` | Nessun valore mancante |
| `rate_ok` | 30–220 bpm |
| `ibi_ok` | Durata entro ±30 % della durata di riferimento |
| `corr_ok` | `template_corr` ≥ **Min. correlation with the reference beat** |
| `amp_ok` | Ampiezza positiva e `amplitude_ratio` tra 1/3 e 3 |
| `clip_ok` | Nessun campione saturato nel battito |
| `valid` | Tutti i controlli superati |
| `variability_raw` | Metrica di variabilità prima della regressione sull'IBI |
| `variability_used` | Metrica usata per H (uguale a `variability_raw` quando la regressione è disattivata) |
| `H` | Variabilità normalizzata (in modalità centrata, il valore letto N_past//2 battiti più avanti) |
| `overlap` | Sovrapposizione O |
| `warmup` | True finché H è tenuto a 0 |
| `local_ibi` | Mediana delle durate degli ultimi 10 battiti validi (s) |
| `T_beats_macro`, `T_beats_time_volume`, `T_beats_derivatives` | Lunghezza della finestra di ogni categoria a questo battito (battiti validi) |

I battiti non validi riportano, nelle colonne della variabilità, i valori dell'ultimo battito valido.

Come leggere la tabella:

* `ibi_ok` falso: di solito un battito ectopico, il battito che lo segue, oppure un picco mancato o in più.
* `corr_ok` falso: un artefatto o un battito di forma atipica.
* `amp_ok` falso: movimento, un cambiamento improvviso della perfusione o del contatto del sensore.
* `clip_ok` falso: saturazione del sensore o dell'amplificatore.
* `rate_ok` falso: quasi sempre un errore di rilevamento.
* `duration_s` dei battiti validi è una serie di intervalli del polso. La variabilità della frequenza del polso
  non coincide con la variabilità della frequenza cardiaca, soprattutto durante stress, movimento e bassa
  perfusione (Schäfer & Vagedes 2013).

---

## 10. Scheda: Variability & windows

Quattro righe condividono l'asse del tempo (tempi dei picchi dei battiti):

1. **Morphological variability (metric)**: `variability_used`, con una linea punteggiata al livello di H = 0 e
   una linea tratteggiata al livello di H = 1.
2. **Normalized variability H**: H, con una linea tratteggiata a θ. Gli intervalli di warm-up sono ombreggiati
   nelle righe 1 e 2.
3. **Window length T_w (valid beats)**: la lunghezza della finestra di ogni categoria (linee a gradini), con
   marcatori ai battiti in cui è stata prodotta una finestra.
4. **Overlap O_w**: la sovrapposizione.

Come leggerlo:

* A riposo, `drift_z` oscilla intorno a 1 e H resta vicino a 0. Le finestre sono a `T_max` e i marcatori sono
  radi.
* Un cambiamento di forma persistente alza `drift_z` per circa 2N battiti, con un massimo circa N battiti dopo
  il cambiamento (trailing). La lunghezza della finestra cala subito e poi torna a crescere di un battito per
  battito (con il limite di crescita di default).
* Gli eventi che possono cambiare la forma dell'impulso includono i cambi di postura, l'inizio e la fine di
  compiti mentali o fisici, gli stimoli vasoattivi (freddo, dolore), i cambiamenti della pressione di contatto
  del sensore e gli artefatti da movimento che hanno superato i controlli di qualità. Confronta i picchi di H
  con i marcatori del protocollo e con gli intervalli scartati.
* Un H periodico alla frequenza respiratoria significa che N_past è troppo breve per lo schema respiratorio
  ([sezione 4](#n_past)).

---

## 11. Scheda: Features

### Controlli

* **Features**: le feature da tracciare, mostrate come "nome [unità]". Default: heart_rate, peak_amplitude,
  crest_time, sdptg_b_a. Sono elencate solo le feature con almeno un valore finito: l'indice di rigidità
  richiede **Subject height**, il PAT richiede un ECG.
* **Show per-beat values**: punti grigi, uno per ogni battito valido.
* **Fixed-window baselines**: le sequenze di finestre fisse da sovrapporre (tutte per default).

### Grafico

Ogni feature è tracciata in un proprio pannello, intitolato "feature [unit] — category: cat".

* **adaptive (category)**: la feature del template di ogni finestra adattiva, a `t_anchor` (il tempo del picco
  del battito di ancoraggio). Con finestre trailing la linea è disegnata a gradini (il valore resta fino alla
  finestra successiva); con finestre centrate i punti sono uniti da segmenti. Passando con il mouse si vedono
  `n_valid` e `T_w`.
* **adaptive CI95**: una banda di ± `1.96 · DS / √n` dei valori per battito della finestra.
* **fixed_Ws**: la feature del template di ogni finestra fissa.
* **ground truth** (solo dati sintetici): la feature senza rumore (linea nera punteggiata).

Sul CI95:

* È l'intervallo al 95 % della **media dei valori per battito**, mentre il valore tracciato è la feature del
  **template**. Per feature non lineari e punti fiduciali fragili le due cose possono differire.
* Usa 1.96 per ogni n. Per finestre piccole il quantile t di Student è maggiore (4.30 con n = 3, 2.26 con
  n = 10), quindi la banda è troppo stretta nelle finestre brevi.
* I battiti successivi non sono indipendenti (modulazione respiratoria e vasomotoria), quindi l'n effettivo è
  minore del numero di battiti.

Considera la banda come un limite inferiore dell'incertezza.

### Adaptive windows table

Una riga per ogni finestra adattiva, tutte le categorie insieme. Colonne:

* `category`, `window_id`;
* `t_anchor`: tempo del picco del battito di ancoraggio;
* `t_start`, `t_end`: onset del primo battito e fine dell'ultimo battito;
* `n_valid`: numero di battiti mediati;
* `template_sqi`: correlazione media dei battiti con il template;
* `anchor_beat`, `T_beats`, `overlap`, `H`;
* ogni feature e il suo `_ci95`. Le colonne delle altre categorie sono vuote (NaN).

> **Suggerimento.** Le finestre brevi sono prodotte molto più spesso di quelle lunghe (una per battito durante
> un transitorio, una ogni ~25 battiti a riposo). Una semplice media o mediana di questa tabella è quindi
> sbilanciata verso i transitori. Per i riassunti nel tempo, usa `features_grid.csv` dalla scheda **Export**,
> che contiene i valori adattivi su una griglia regolare a 1 Hz.

### Per-beat features table

Una riga per ogni battito rilevato (compresi i battiti non validi; filtra su `valid`): `beat_index`, `t_peak`,
`valid`, e le 33 feature misurate su quel singolo battito con la stessa procedura dei template.

> **Suggerimento.** I valori del template e quelli per battito sono stimatori diversi. Sui dati di esempio,
> c/a, d/a, e/a e l'indice di invecchiamento per battito differiscono sistematicamente dai valori del template,
> e le feature per battito legate all'incisura di BIDMC 01 sono instabili
> ([sezione 19](#19-valori-di-riferimento-dei-dati-di-esempio)). Confronta valori del template con valori del
> template, e valori per battito con valori per battito.

---

## 12. Scheda: Templates

**Category** sceglie la categoria. **Templates shown** (2–30, default 12) sceglie quanti template disegnare,
distribuiti uniformemente nel tempo e colorati dal blu (iniziali) al rosso (finali). La legenda riporta il
tempo dell'ancora e il numero di battiti ("t=XXs (n=NN)").

Come si costruisce un template:

1. I battiti della finestra sono campionati su una griglia a 500 Hz per interpolazione B-spline cubica del
   segnale di morfologia.
2. Sono allineati sul punto di massima pendenza sistolica (precisione inferiore al campione).
3. Si rimuove la linea di base retta da piede a piede di ogni battito, e ogni battito è diviso per la sua
   ampiezza.
4. I battiti sono mediati, e la media è moltiplicata per l'ampiezza media.

Il template conserva la forma media e l'ampiezza media, sulla scala dei tempi reale (i battiti non sono stirati
a una durata comune). La sua durata è il tempo mediano dal piede al punto di massima pendenza più il tempo
mediano da quel punto alla fine del battito. L'asse y è l'ampiezza in u.a. sopra la linea di base da piede a
piede.

> **Problema noto.** L'asse del tempo ("Time from foot (s)") è dilatato di un fattore 500 / frequenza di
> elaborazione (4 volte a 125 Hz). I valori delle feature non ne sono influenzati. Leggi i tempi dalle
> feature, non da questo asse ([sezione 20](#20-problemi-noti)).

Cosa osservare:

* Un picco sistolico, poi un'incisura dicrota oppure un'inflessione tardo-sistolica o diastolica.
* Con l'età e l'irrigidimento arterioso, l'onda diastolica tende a fondersi con la parte sistolica e
  l'incisura diventa meno visibile (Millasseau et al. 2002, 2006). In questo caso il programma usa il punto di
  inflessione.
* Template che cambiano gradualmente dal blu al rosso indicano una lenta deriva della forma. Un salto
  improvviso indica un transitorio o un artefatto.

---

## 13. Scheda: Benchmark

Questa scheda funziona solo con **Synthetic demo**. Con le altre sorgenti mostra una nota che suggerisce i
comandi `adaptive-ppg benchmark` e `adaptive-ppg sweep` per studi più ampi.

Il riferimento vero (ground truth) è lo stesso scenario senza rumore, deriva, respirazione, jitter o artefatti,
analizzato alla frequenza nativa. Le sue feature per battito sono interpolate linearmente sulla griglia a 1 Hz
dell'analisi. Ogni stimatore è confrontato con esso dopo `max(60 s, la finestra fissa più lunga)`. Stimatori:

* `adaptive`: le finestre adattive;
* `fixed_Ws`: le finestre fisse;
* `per_beat`: l'ultimo battito valido (mantenimento di ordine zero);
* `ema_10s`: una media mobile esponenziale causale dei valori per battito, con costante di tempo di 10 s.

Le feature valutate sono quelle il cui riferimento vero varia abbastanza nello scenario di stress (DS almeno
pari al 2 % del valore assoluto medio) e che sono disponibili almeno per l'80 % del tempo. Le scale sono le DS
del riferimento vero nello scenario di stress (seed 0).

### Tabella riassuntiva

| Colonna | Significato |
| :--- | :--- |
| `estimator` | Nome dello stimatore |
| `median_nrmse` | Mediana sulle feature di RMSE / DS del riferimento vero (più basso è meglio) |
| `median_stationary_nsd` | Mediana sulle feature della DS della stima entro i segmenti stazionari, normalizzata con la stessa scala (rumore a riposo; più basso è meglio) |
| `median_latency_s` | Mediana sulle feature del ritardo (s) dell'attraversamento del 50 % alle transizioni |
| `mean_coverage` | Quota media della griglia in cui lo stimatore ha un valore |

La tabella è ordinata per `median_nrmse`. I segmenti stazionari dello scenario di stress sono 60–120 s e
140–240 s. Il recupero non è incluso perché non si stabilizza entro 120 s (costante di tempo 25 s). La latenza
di una transizione è contata solo per le feature che cambiano di almeno 0.5 DS in quella transizione. I suoi
livelli prima e dopo sono le medie dei 30 s prima dell'evento e dei 30 s prima del confine successivo.

**Metric** (nrmse, stationary_nsd, latency_s, bias) disegna un box plot per ogni stimatore, sui valori delle
feature, con gli stimatori ordinati per mediana. Il bias è nelle unità della feature, quindi non si può
confrontare tra feature.

**Per-feature metrics** elenca `coverage`, `rmse`, `bias`, `nrmse`, `stationary_sd`, `stationary_nsd`,
`latency_s` e `n_transitions` per ogni stimatore e feature.

Esempio (stress, seed 1, 125 Hz, default):

| Stimatore | median_nrmse | median_stationary_nsd | median_latency_s |
| :--- | :---: | :---: | :---: |
| fixed_10s | 0.340 | 0.127 | 8.0 |
| ema_10s | 0.361 | 0.201 | 10.0 |
| adaptive | 0.422 | 0.090 | 8.5 |
| per_beat | 0.494 | 0.474 | −2.0 |
| fixed_30s | 0.619 | 0.342 | 23.0 |

La copertura è 1.0 per tutti gli stimatori. Le finestre adattive hanno il rumore più basso a riposo e una
latenza simile a quella delle finestre fisse di 10 s. In questo scenario una finestra di 10 s ha un errore
complessivo più basso. La latenza di −2 s del valore per battito non è un vero anticipo: è prodotta dal rumore
all'attraversamento.

> **Avvertenze.** Il modello sintetico è semplice. Il benchmark confronta gli stimatori in condizioni note; non
> valida le feature sulla fisiologia reale. Il README riporta i risultati su tre seed e per lo scenario di
> riposo.

---

## 14. Scheda: Export

* **Include ensemble templates**: aggiunge alle esportazioni i template di ogni categoria.
* **Download ZIP (CSV + metadata.json)**:

| File | Contenuto |
| :--- | :--- |
| `features_grid.csv` | Feature adattive su una griglia a 1 Hz, dal primo secondo intero dopo il primo battito all'ultimo battito. Trailing: il valore della finestra più recente (vuoto prima della prima). Centrata: interpolazione lineare tra le ancore. |
| `windows.csv` | La tabella delle finestre adattive ([sezione 11](#adaptive-windows-table)) |
| `beats.csv` | La tabella della qualità dei battiti ([sezione 9](#beat-quality-table)) |
| `features_per_beat.csv` | Le feature per battito |
| `feature_dictionary.csv` | Nome, categoria e unità di ogni feature |
| `windows_fixed_<W>s.csv` | Una tabella per ogni finestra fissa |
| `templates_<category>.csv` | I template (solo con **Include ensemble templates**) |
| `metadata.json` | Metadati dell'esecuzione (vedi sotto) |

* **Download Excel workbook**: le stesse tabelle come fogli (nomi accorciati a 31 caratteri) più un foglio
  `metadata`. Richiede il pacchetto `openpyxl`; altrimenti compare "Excel export needs openpyxl.".
* **Download configuration (JSON)**: `config.json`, con tutti i valori dei parametri. Si può riusare con
  `adaptive-ppg run --config config.json`. **Non** contiene il nome del preset né l'elenco delle finestre
  fisse (contiene `fixed_overlap_frac` e un `fixed_window_sec` di 30 che la GUI non usa).
* **Run metadata**: mostra `metadata.json`:
  * pacchetto, versione, revisione git (con `-dirty` quando i file tracciati hanno modifiche non committate),
    ora di creazione (UTC), rilevatore dei picchi;
  * input: canale, formato, frequenza di campionamento, durata, e i dettagli del caricamento (colonna
    risolta, separatore, separatore decimale, colonna del tempo, frequenza stimata, deviazione massima dei
    tempi, origine della frequenza, avviso sulla frequenza);
  * summary: numero di battiti e di battiti validi, frequenza di elaborazione, polarità e asimmetria, frequenza
    di taglio effettiva del passa-basso, numero di finestre per categoria, pendenza della regressione sull'IBI;
  * avvisi e configurazione completa.

---

## 15. Guida alle feature

Tutte le feature sono misurate su un template (o su un singolo battito) campionato a 500 Hz. Il template è
smussato con un filtro di Savitzky–Golay (70 ms, ordine 3; Savitzky & Golay 1964), che fornisce anche la
derivata prima (`dy`) e la derivata seconda (`ddy`). Notazione: onset `o`, picco sistolico `p` (massimo tra
onset e fine, con raffinamento parabolico), fine `e`, linea di base `base` = valore smussato all'onset,
`amp` = valore al picco − `base`. I tempi usano posizioni con precisione inferiore al campione.

Le letture fisiologiche qui sotto riassumono la letteratura citata. Descrivono associazioni riportate in
popolazioni e configurazioni specifiche, non misure calibrate.

### Macro (ampiezza e forma dell'impulso)

| Feature | Definizione come implementata | Unità | Lettura fisiologica | Avvertenze |
| :--- | :--- | :---: | :--- | :--- |
| `peak_amplitude` | `amp` | u.a. | Variazione pulsatile (AC) del volume ematico nel tessuto sotto il sensore. Diminuisce con la vasocostrizione locale (per esempio attivazione simpatica, freddo, dolore) e aumenta con la vasodilatazione (Allen 2007). | Dipende da sede, lunghezza d'onda, pressione di contatto e guadagno del dispositivo. Non confrontabile tra soggetti o dispositivi. Misurata sul segnale filtrato. |
| `pulse_width_10`, `_25`, `_50`, `_75` | Tempo tra l'attraversamento in salita e quello in discesa del 10, 25, 50, 75 % di `amp` | s | È stato riportato che la larghezza a metà altezza è correlata con la resistenza vascolare sistemica (Awad et al. 2007). | Si accorcia con la frequenza cardiaca. La larghezza al 10 % risente dell'onda diastolica e dell'incisura. |
| `total_area` | Area dell'impulso sopra `base`, dall'onset alla fine | u.a.·s | Scala con il volume pulsatile e con la durata del battito | Stessa dipendenza dall'ampiezza di `peak_amplitude` |
| `systolic_area`, `diastolic_area` | Area dall'onset al picco sistolico, e dal picco sistolico alla fine | u.a.·s | — | Divise al **picco** sistolico, non all'incisura |
| `area_ratio` | `diastolic_area / systolic_area` | – | Equilibrio tra la parte tardiva e quella iniziale dell'impulso | Dipende fortemente dal duty cycle, quindi dalla frequenza cardiaca (quando la frequenza cardiaca aumenta, la diastole si accorcia più della sistole) |
| `inflection_point_area_ratio` (IPA) | Area dall'incisura alla fine divisa per l'area dall'onset all'incisura | – | Proposto come indice legato alla resistenza periferica totale (Wang et al. 2009) | Dipende dalla posizione dell'incisura. Rumoroso: nel benchmark del README il suo NRMSE supera 1 con tutti gli stimatori. |
| `inflection_amplitude` | Valore smussato al punto diastolico − `base` | u.a. | Altezza dell'onda diastolica (riflessa) | Il punto diastolico è un picco quando esiste un'incisura, altrimenti un'inflessione |
| `reflection_index` (RI) | `inflection_amplitude / amp` | – | Riflette soprattutto il tono delle piccole arterie. Diminuisce con vasodilatatori come la nitroglicerina (gliceril trinitrato) e con il β2-agonista salbutamolo (Chowienczyk et al. 1999; Millasseau et al. 2006). | Rapporto, non percentuale (moltiplica per 100 per confrontarlo con studi in %). Instabile per battito quando non c'è un'incisura chiara. |
| `max_systolic_slope` | Massimo di `dy` tra onset e picco | u.a./s | Ripidità della salita | Scala con l'ampiezza |
| `max_decay_slope` | Minimo di `dy` tra picco e fine (negativo) | u.a./s | Ripidità della discesa | Scala con l'ampiezza |
| `slope_ratio` | abs(`max_systolic_slope / max_decay_slope`) | – | Asimmetria tra salita e discesa | — |
| `dc_level` | Media del segnale grezzo sul battito (ricampionato, non filtrato, non invertito) | u.a. | Con l'intensità luminosa grezza: il livello di luce non pulsatile (tessuto, sangue venoso e sangue arterioso non pulsatile, e intensità della sorgente) | Privo di significato quando il dispositivo rimuove o riscala la componente DC |
| `perfusion_index` (PI) | `100 · amp / dc_level`; nelle finestre, ampiezza del template / `dc_level` medio | % | Rapporto tra segnale pulsatile e non pulsatile. Valori bassi indicano una scarsa perfusione periferica (Lima et al. 2002). | Ha significato solo per un'intensità grezza accoppiata in DC, dallo stesso canale. Qui la componente AC è misurata sul segnale filtrato, quindi i valori non coincidono con il PI dei pulsossimetri. Non fisiologico per i dati sintetici e per BIDMC. |

### Time_volume (tempi)

| Feature | Definizione come implementata | Unità | Lettura fisiologica | Avvertenze |
| :--- | :--- | :---: | :--- | :--- |
| `pulse_duration` | Dall'onset alla fine (onset successivo) | s | Intervallo del polso | Per un template, l'intervallo mediano dal piede alla massima pendenza più l'intervallo mediano dalla massima pendenza alla fine dei suoi battiti, vicino alla durata mediana dei battiti |
| `heart_rate` | 60 / `pulse_duration` | bpm | Frequenza del polso | Frequenza del polso, non frequenza cardiaca da ECG (Schäfer & Vagedes 2013) |
| `crest_time` | Dall'onset al picco sistolico | s | Durata della salita. Usato con l'indice di rigidità per classificare la rigidità arteriosa (Alty et al. 2007). | Dipende da **Onset method** (più breve con `tangent`) |
| `decay_time` | Dal picco sistolico alla fine | s | Comprende la diastole, che è la parte che si accorcia di più quando la frequenza cardiaca aumenta | — |
| `duty_cycle` | `crest_time / pulse_duration` | – | Frazione del ciclo occupata dalla salita | Aumenta con la frequenza cardiaca |
| `notch_time` | Dall'onset all'incisura dicrota; senza incisura, dall'onset al massimo di `ddy` tra il picco e il punto di inflessione | s | Riferimento temporale della tarda sistole nel sito di misura | Senza validazione, non è una misura del tempo di eiezione ventricolare sinistra |
| `delta_t_dvp` (ΔT_DVP) | Dal picco sistolico al picco diastolico (o al punto di inflessione) | s | Tempo tra l'onda diretta e quella riflessa. Si accorcia con la rigidità delle grandi arterie e con l'età (Millasseau et al. 2002). | Richiede un punto diastolico risolvibile |
| `stiffness_index` (SI) | Altezza del soggetto / `delta_t_dvp` | m/s | Indice di rigidità delle grandi arterie, correlato con la velocità dell'onda di polso carotido-femorale (Millasseau et al. 2002) | Richiede **Subject height**. In assenza di picco diastolico si usa il punto di inflessione, come in Millasseau et al. (2002). |
| `pat_foot`, `pat_peak` | Dal picco R dell'ECG precedente all'onset del PPG, e al picco sistolico | s | Tempo di arrivo dell'impulso = periodo di pre-eiezione + tempo di transito dell'impulso (Mukkamala et al. 2015). Si accorcia con una pressione arteriosa più alta (tempo di transito più breve) e con una contrattilità maggiore (periodo di pre-eiezione più breve). | Richiede un ECG. Misurato su ogni battito alla frequenza di elaborazione; i valori delle finestre sono medie dei valori per battito. I valori fuori da 0.05–0.6 s (piede) sono scartati. I ritardi del dispositivo tra il canale ECG e quello pletismografico aggiungono un offset. |

### Derivatives (derivata seconda, SDPTG)

Le onde della SDPTG (Takazawa et al. 1998; Elgendi 2012) sono trovate su `ddy` come segue:

* **a**: massimo di `ddy` da 50 ms prima dell'onset al punto di massima pendenza;
* **b**: minimo di `ddy` dopo a, fino al picco sistolico;
* **e**: massimo di `ddy` dopo il picco, entro il primo 60 % del battito;
* **c**: il massimo locale più alto tra b ed e, con una prominenza di almeno il 5 % di a;
* **d**: il minimo locale più basso tra c ed e, con la stessa prominenza. Se d non viene trovato, si scarta
  anche c.

I rapporti sono calcolati solo quando a > 0 e b < 0. Quando c o d non si possono risolvere (cosa frequente nei
soggetti anziani), c/a, d/a e l'indice di invecchiamento sono vuoti (NaN).

| Feature | Definizione | Lettura fisiologica |
| :--- | :--- | :--- |
| `sdptg_b_a` | b/a | Negativo. Aumenta (verso 0) con l'età e la rigidità arteriosa (Takazawa et al. 1998). |
| `sdptg_c_a` | c/a | Diminuisce con l'età (Takazawa et al. 1998) |
| `sdptg_d_a` | d/a | Diminuisce con l'età e cambia con gli agenti vasoattivi (Takazawa et al. 1998) |
| `sdptg_e_a` | e/a | Diminuisce con l'età (Takazawa et al. 1998) |
| `aging_index` | (b − c − d − e)/a | Indice di invecchiamento della SDPTG; aumenta con l'età (Takazawa et al. 1998) |
| `aging_index_be` | (b − e)/a | Indice alternativo quando c e d non si possono risolvere (rassegna in Elgendi 2012) |

> **Importante.** I rapporti della SDPTG dipendono fortemente dallo smussamento della derivata, dalla
> frequenza di taglio del passa-basso e dalla frequenza di campionamento. Sullo stesso segnale sintetico il
> README riporta b/a = −0.65 con una finestra di smussamento di 50 ms, −1.12 con 70 ms e −1.64 con 80 ms.
> Confronta i valori della SDPTG solo tra analisi con impostazioni identiche, e confrontali con le norme
> pubblicate solo se l'elaborazione è la stessa. Per una rassegna di questi indici e dei loro limiti, vedi
> Charlton et al. (2022).

---

## 16. Impostazioni per caso d'uso

Questi sono punti di partenza. Verifica i risultati con gli strumenti della
[sezione 17](#17-lista-dei-controlli).

### Registrazioni a riposo e caratterizzazione della linea di base

* **Preset** `rest`. Scrivi a mano **T_max (s)** 60 / 60 / 90 in **3. Windows per category**
  ([sezione 20](#20-problemi-noti)).
* **H normalization** fixed con `drift_z`.
* **Window anchor** `centered` per l'analisi offline.
* **Fixed windows** `30, 60` per confronto.
* Riassumi ogni registrazione con la mediana di `features_grid.csv`.
* Con il preset di default, le registrazioni di esempio a riposo hanno H > θ nel 21–26 % dei battiti
  (`Sample1`, `Sample2`). Il preset `rest` (θ = 0.5) fa ridurre le finestre solo per cambiamenti più grandi.

### Protocolli autonomici acuti e di stress

Esempi: aritmetica mentale, cold pressor test, manovra di Valsalva, handgrip, ortostatismo attivo o tilt test.

* **Preset** `acute_stress`. Scrivi a mano **T_max (s)** 20 / 20 / 40.
* H non vede i cambiamenti di sola ampiezza o di soli tempi. Per seguire l'ampiezza dell'impulso, l'indice di
  perfusione o la frequenza del polso entro pochi secondi, abbassa ancora il **T_max (s)** di macro e
  time_volume (per esempio 10–15 s).
* Evita la normalizzazione **calibration**, a meno che i primi 60 s non siano una linea di base pulita.
  Preferisci **fixed** con `drift_z`.
* Un **N_past** più piccolo (per esempio 6–8) reagisce più in fretta, ma è meno sensibile e più influenzato dal
  respiro.
* Per medie sincronizzate agli eventi, usa `centered` offline, ricordando che H ritarda comunque di circa N/2
  battiti. Usa `trailing` quando conta la causalità, e ricorda il ritardo di circa T/2.
* Aggiungi un ECG per il PAT. La sua componente di pre-eiezione si accorcia con l'attivazione simpatica, quindi
  i cambiamenti del PAT combinano effetti cardiaci e vascolari.
* Aggiungi **Fixed windows** `10` per confronto.

### Analisi causale o simile al tempo reale

* **Window anchor** `trailing`. **H normalization** fixed, calibration o rolling (non global).
* La regressione sull'IBI è causale. Solo per i primi 5 battiti il riferimento di qualità usa i primi 30
  battiti.
* La dashboard analizza file interi. Per l'uso in tempo reale, il pacchetto fornisce
  `StreamingAdaptiveEngine`, che riceve un battito alla volta e restituisce le finestre adattive man mano che
  si chiudono (vedi il README).

### Dispositivi indossabili da polso e a bassa frequenza (circa 64 Hz)

* **Preset** `wearable_64hz`. Scrivi a mano 15 nel **T_crit** di derivatives.
* Lascia **Processing rate (Hz)** a 125. Il sovracampionamento non aggiunge banda: le frequenze di taglio sono
  comunque limitate al 45 % di 64 Hz (28.8 Hz).
* Sotto 100 Hz le feature della derivata non sono consigliate (README: differenza di 0.3 DS a 64 Hz). Sotto
  50 Hz compare un avviso.
* La forma dell'impulso dipende dal sito di misura. Le feature basate sull'incisura e i valori di riferimento
  degli studi al dito non si trasferiscono direttamente al polso. Gli artefatti da movimento sono più
  frequenti: controlla la quota di battiti validi.

### Respirazione lenta o guidata, aritmia sinusale respiratoria, biofeedback

* **N_past** ≥ 2 × frequenza cardiaca / frequenza respiratoria: per 6 atti/min a 60 bpm, N ≥ 20.
* Valuta **Regress out IBI-driven variability** se H segue il respiro.
* Le finestre brevi mostrano la modulazione respiratoria delle feature. Le finestre lunghe la eliminano con la
  media.

### Battiti ectopici, aritmie e fibrillazione atriale

* I controlli di qualità (`ibi_ok`, `corr_ok`) di solito scartano i battiti ectopici e i battiti che li
  seguono. Non abbassare **Min. correlation with the reference beat** per conservarli.
* Nella fibrillazione atriale gli intervalli sono irregolari e l'ampiezza dell'impulso varia con l'intervallo
  precedente. Molti battiti sono scartati e la media d'insieme non ha senso. Il metodo non è pensato per questo
  ritmo: analizza solo i segmenti in ritmo sinusale.

### Esportazioni da monitor clinici (MIMIC, BIDMC)

* La curva pletismografica è già filtrata, scalata e mostrata con il volume verso l'alto, di solito a 125 Hz.
  I default funzionano.
* `perfusion_index` e `dc_level` non hanno significato (BIDMC 01 dà un PI di circa 37 %).
* Il PAT può includere ritardi di elaborazione tra i canali ECG e pletismografico del monitor. Il valore di
  BIDMC 01 (`pat_foot` ≈ 0.53 s) è più lungo dei valori di solito riportati per il dito, circa 0.2–0.4 s.
  Interpreta i cambiamenti all'interno di una registrazione piuttosto che i valori assoluti.

### Bradicardia e tachicardia

* Sotto circa 45 bpm, abbassa **Morphology high-pass (Hz)** a 0.3 Hz. I battiti sotto 30 bpm sono sempre
  scartati.
* Il rilevatore di Elgendi non trova picchi più vicini di 0.3 s, quindi è limitato a 200 bpm. I battiti sopra
  220 bpm sono sempre scartati. I default non sono stati validati per neonati o bambini piccoli.
* **T_max (s)** è in secondi: a frequenze cardiache alte la stessa finestra contiene più battiti.

### Confronto tra soggetti, gruppi o sessioni

* Usa la stessa configurazione (salva il JSON) e le stesse finestre fisse.
* Mantieni uguali la frequenza di elaborazione, il passa-basso, il metodo di onset e lo smussamento della
  derivata.
* Confronta valori del template con valori del template, e valori per battito con valori per battito.
* Riassumi con valori pesati nel tempo (`features_grid.csv`), non con la tabella delle finestre.
* Ampiezza, aree, pendenze, DC e PI sono in unità arbitrarie e dipendono dal dispositivo e dal contatto.
  Preferisci i cambiamenti entro il soggetto, oppure le feature adimensionali e temporali.
* Le feature temporali e il rapporto tra aree dipendono dalla frequenza cardiaca. Valuta di usare la frequenza
  cardiaca come covariata.
* Riporta la quota di battiti validi, e stabilisci in anticipo una regola di esclusione.

---

## 17. Lista dei controlli

1. La riga informativa mostra il canale, la frequenza di campionamento e la durata giusti. Nessun avviso di
   discrepanza della frequenza.
2. L'impulso filtrato punta verso l'alto. Altrimenti, imposta a mano **Signal polarity**.
3. Picchi e onset sono sui punti giusti quando ingrandisci alcuni battiti.
4. **Valid beats** è sopra circa il 90 %. Gli intervalli scartati corrispondono ad artefatti visibili.
5. **Median template SQI** è vicino a 1 (i dati di esempio danno 0.985–0.993 per macro).
6. Il warm-up è breve rispetto alla registrazione (19 battiti validi con i default).
7. **Beats with H > θ** è plausibile: basso a riposo, più alto durante i compiti.
8. I template hanno una forma plausibile. Controlla se l'incisura è visibile.
9. Per le feature SDPTG: controlla quanti battiti hanno c/a e d/a vuoti in **Per-beat features**.
10. Per il PAT: i valori non si accumulano ai limiti (0.05 o 0.6 s), e il rilevamento dei picchi R è
    plausibile.

---

## 18. Cosa riportare

* Versione del pacchetto e revisione git (`metadata.json`).
* Il preset e ogni parametro modificato (`config.json`), più le finestre fisse, che non sono nel JSON.
* Frequenza di elaborazione, bande dei filtri (e frequenza di taglio effettiva del passa-basso), metodo di
  onset, rilevatore dei picchi.
* Smussamento della derivata (Savitzky–Golay 70 ms, ordine 3, salvo modifiche) e griglia dei template a
  500 Hz.
* Soglie di qualità e quota di battiti validi.
* Metrica di variabilità, normalizzazione di H, θ, pendenza, N_past, ancoraggio e limiti delle finestre per
  categoria.
* Se i valori riportati sono valori del template o per battito, e come sono stati riassunti nel tempo.

---

## 19. Valori di riferimento dei dati di esempio

Impostazioni di default, elaborazione a 125 Hz, finestre fisse di 10 e 30 s. "Template" è la mediana sulle
finestre adattive della categoria. "Per battito" è la mediana sui battiti validi. Questi valori descrivono le
tre registrazioni di esempio. Non sono valori normativi.

| | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01_Signals.csv` |
| :--- | :---: | :---: | :---: |
| Battiti validi | 335 / 335 | 276 / 277 | 705 / 718 |
| Frequenza del polso (bpm) | 71.1 | 73.7 | 90.9 |
| `drift_z` mediano | 1.18 | 1.06 | 0.83 |
| Battiti con H > θ | 25.9 % | 21.0 % | 7.2 % |
| Finestra macro, mediana su tutti i battiti (p10–p90 dopo il warm-up) | 18 (4–33) | 17 (5–34) | 36 (13–44) |
| Finestra macro, mediana sulle finestre | 8 | 10.5 | 29 |
| Numero di finestre macro | 47 | 32 | 40 |
| Finestra derivatives, mediana su tutti i battiti / sulle finestre | 25 / 16.5 | 25 / 24 | 60 / 46 |
| SQI del template, macro / derivatives | 0.985 / 0.958 | 0.989 / 0.982 | 0.993 / 0.993 |
| Warm-up | 19 battiti (15.8 s) | 20 battiti (15.0 s) | 19 battiti (11.9 s) |

Le mediane su tutti i battiti includono il warm-up, come il KPI **Median macro window**. Il README riporta le
mediane dopo il warm-up (macro: 17, 18 e 35 battiti; derivatives: 24, 26 e 59 battiti).

Feature (template / per battito):

| Feature | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01_Signals.csv` |
| :--- | :---: | :---: | :---: |
| `crest_time` (s) | 0.095 / 0.096 | 0.093 / 0.095 | 0.111 / 0.111 |
| `duty_cycle` | 0.112 / 0.115 | 0.116 / 0.116 | 0.167 / 0.167 |
| `notch_time` (s) | 0.285 / 0.302 | 0.261 / 0.294 | 0.224 / 0.224 |
| `delta_t_dvp` (s) | 0.242 / 0.247 | 0.249 / 0.247 | 0.155 / 0.339 |
| `reflection_index` | 0.657 / 0.665 | 0.628 / 0.665 | 0.289 / −0.033 |
| `pulse_width_50` (s) | 0.429 / 0.424 | 0.404 / 0.416 | 0.166 / 0.166 |
| `area_ratio` | 5.63 / 5.29 | 5.36 / 5.31 | 1.60 / 1.64 |
| `inflection_point_area_ratio` | 0.81 / 0.70 | 0.88 / 0.71 | 0.147 / 0.161 |
| `slope_ratio` | 6.65 / 4.99 | 6.38 / 4.78 | 1.82 / 1.68 |
| `sdptg_b_a` | −0.724 / −0.732 | −0.702 / −0.709 | −0.874 / −0.898 |
| `sdptg_c_a` | 0.058 / 0.116 | 0.080 / 0.119 | −0.413 / −0.307 |
| `sdptg_d_a` | −0.099 / −0.156 | −0.105 / −0.165 | −0.484 / −0.527 |
| `sdptg_e_a` | 0.100 / 0.190 | 0.093 / 0.187 | 0.451 / 0.499 |
| `aging_index` | −0.787 / −0.876 | −0.777 / −0.860 | −0.443 / −0.580 |
| `aging_index_be` | −0.831 / −0.922 | −0.789 / −0.904 | −1.326 / −1.397 |
| `perfusion_index` (%) | 0.337 / 0.381 | 0.353 / 0.357 | 36.7 / 37.6 (senza significato) |

BIDMC 01 con la derivazione ECG II: `pat_foot` ≈ 0.525 s, `pat_peak` ≈ 0.637 s (mediane).

Osservazioni:

* c/a, d/a e l'indice di invecchiamento per battito sono vuoti nel 15 % (`Sample1`), nel 16 % (`Sample2`) e
  nel 22 % (BIDMC 01) dei battiti.
* c/a, d/a, e/a e l'indice di invecchiamento per battito differiscono sistematicamente dai valori del
  template: la media smussa le piccole onde c e d.
* In BIDMC 01 ΔT_DVP e RI per battito sono instabili (l'RI mediano per battito è negativo). Usa i valori del
  template.
* Tassi di superamento dei singoli controlli in BIDMC 01: `rate_ok` 99.6 %, `ibi_ok` 98.5 %, `corr_ok`
  99.3 %, `amp_ok` 99.0 %.
* I valori molto diversi di `area_ratio` e `pulse_width_50` di BIDMC 01 riflettono un diverso contorno
  dell'impulso (frequenza cardiaca più alta, paziente critico, filtraggio del monitor). Non vanno letti come
  una differenza tra sensori.

---

## 20. Problemi noti

Questi problemi sono stati trovati verificando questa guida rispetto al codice della versione 0.2.0. Non sono
corretti in quella versione.

1. **I preset non aggiornano la sezione 3 della barra laterale.** I controlli **T_min (beats)**,
   **T_max (s)** e **T_crit** mantengono i valori precedenti quando cambi **Preset** dopo la prima esecuzione.
   Scrivi a mano i valori di `rest`, `acute_stress` e `wearable_64hz` per **3. Windows per category**
   ([sezione 4](#preset)). Gli altri controlli seguono il preset, e i valori senza controllo (limite di
   crescita, smussamento della derivata) si applicano sempre.
2. **Asse del tempo dei template.** Nella scheda **Templates**, e nella colonna `t_from_foot_s` dei template
   esportati, l'asse del tempo è diviso per la frequenza di elaborazione invece che per la frequenza dei
   template, 500 Hz. È dilatato di un fattore 500 / frequenza di elaborazione (4 volte a 125 Hz; corretto
   quando si elabora a 500 Hz). I valori delle feature non ne sono influenzati.
3. **README: normalizzazione dei battiti.** Il README dice che ogni battito è scalato in [0, 1]. Il codice
   sottrae la retta dall'onset alla fine e scala il picco a 1, quindi sono possibili valori sotto 0.
4. **README: finestre centrate.** Il README dice che `centered` dà finestre a ritardo zero. Le finestre delle
   feature sono simmetriche, ma H reagisce comunque circa N_past/2 battiti dopo un cambiamento.
5. **Impostazioni non salvate.** `config.json` non contiene il nome del preset né l'elenco delle finestre
   fisse. La scheda **Benchmark** usa sempre il rilevatore di Elgendi, qualunque cosa indichi
   **Peak detector**.

---

## 21. Bibliografia

* Allen J. Photoplethysmography and its application in clinical physiological measurement. *Physiol Meas*
  2007;28:R1–R39.
* Alty SR, et al. Predicting arterial stiffness from the digital volume pulse waveform. *IEEE Trans Biomed
  Eng* 2007;54:2268–2275.
* Awad AA, et al. The relationship between the photoplethysmographic waveform and systemic vascular
  resistance. *J Clin Monit Comput* 2007;21:365–372.
* Charlton PH, et al. Assessing hemodynamics from the photoplethysmogram to gain insights into vascular age:
  a review from VascAgeNet. *Am J Physiol Heart Circ Physiol* 2022;322:H493–H522.
* Chowienczyk PJ, et al. Photoplethysmographic assessment of pulse wave reflection: blunted response to
  endothelium-dependent beta2-adrenergic vasodilation in type II diabetes mellitus. *J Am Coll Cardiol*
  1999;34:2007–2014.
* Elgendi M. On the analysis of fingertip photoplethysmogram signals. *Curr Cardiol Rev* 2012;8:14–25.
* Elgendi M, et al. Systolic peak detection in acceleration photoplethysmograms measured from emergency
  responders in tropical conditions. *PLoS ONE* 2013;8:e76585.
* Goldberger AL, et al. PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for
  complex physiologic signals. *Circulation* 2000;101:e215–e220.
* Julien C. The enigma of Mayer waves: facts and models. *Cardiovasc Res* 2006;70:12–21.
* Li Q, Clifford GD. Dynamic time warping and machine learning for signal quality assessment of pulsatile
  signals. *Physiol Meas* 2012;33:1491–1501.
* Lima AP, Beelen P, Bakker J. Use of a peripheral perfusion index derived from the pulse oximetry signal as a
  noninvasive indicator of perfusion. *Crit Care Med* 2002;30:1210–1213.
* Mehrgardt P, Khushi M, Poon S, Withana A. Pulse Transit Time PPG Dataset (version 1.1.0). PhysioNet 2022.
  doi:10.13026/jpan-6n92.
* Millasseau SC, et al. Determination of age-related increases in large artery stiffness by digital pulse
  contour analysis. *Clin Sci* 2002;103:371–377.
* Millasseau SC, et al. Contour analysis of the photoplethysmographic pulse measured at the finger.
  *J Hypertens* 2006;24:1449–1456.
* Mukkamala R, et al. Toward ubiquitous blood pressure monitoring via pulse transit time: theory and
  practice. *IEEE Trans Biomed Eng* 2015;62:1879–1901.
* Orphanidou C, et al. Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and
  applications to wireless monitoring. *IEEE J Biomed Health Inform* 2015;19:832–838.
* Pimentel MAF, et al. Toward a robust estimation of respiratory rate from pulse oximeters. *IEEE Trans Biomed
  Eng* 2017;64:1914–1923.
* Savitzky A, Golay MJE. Smoothing and differentiation of data by simplified least squares procedures.
  *Anal Chem* 1964;36:1627–1639.
* Schäfer A, Vagedes J. How accurate is pulse rate variability as an estimate of heart rate variability? A
  review on studies comparing photoplethysmographic technology with an electrocardiogram. *Int J Cardiol*
  2013;166:15–29.
* Takazawa K, et al. Assessment of vasoactive agents and vascular aging by the second derivative of
  photoplethysmogram waveform. *Hypertension* 1998;32:365–370.
* Wang L, et al. Noninvasive cardiac output estimation using a novel photoplethysmogram index. *Proc IEEE Eng
  Med Biol Soc (EMBC)* 2009:1746–1749.
