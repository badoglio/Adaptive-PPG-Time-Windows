# Guida all'installazione, passo per passo

[English](INSTALL_en.md) · [Italiano](INSTALL_it.md) · [日本語](INSTALL_ja.md)

Questa guida è pensata per chi non ha mai installato un programma Python. Richiede circa 15 minuti, una
connessione a internet e circa 1 GB di spazio libero su disco.

Alla fine avrai:

* la **dashboard**, una pagina web che gira sul tuo computer e analizza registrazioni PPG;
* il **comando `adaptive-ppg`**, per analizzare i file dal terminale.

I comandi da digitare sono nei riquadri grigi. Scrivili (o copiali e incollali) nel terminale e premi
**Invio**. Dove la guida cambia a seconda del sistema operativo, segui solo la parte che ti riguarda.

---

## Passo 1. Installare Python

Il programma richiede **Python 3.10 o successivo**. Consigliamo **Python 3.12**, la versione su cui è
testato.

### Windows

1. Apri <https://www.python.org/downloads/windows/> e scarica il **Windows installer (64-bit)** di
   Python 3.12.
2. Avvia il file scaricato.
3. **Importante:** nella prima schermata spunta **"Add python.exe to PATH"** in basso, poi fai clic su
   **"Install Now"**.
4. Quando compare "Setup was successful", fai clic su **Close**.

### macOS

1. Apri <https://www.python.org/downloads/macos/> e scarica il **macOS 64-bit universal2 installer** di
   Python 3.12.
2. Apri il file `.pkg` e segui l'installazione.
3. Su macOS il comando si chiama **`python3`**, non `python`. Dove questa guida scrive `python`, digita
   `python3` (solo fino al Passo 5: dentro l'ambiente virtuale funziona anche `python`).

### Linux (Ubuntu / Debian)

Di solito Python è già installato. Installa anche il modulo per gli ambienti virtuali:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip
```

Come su macOS, digita `python3` al posto di `python` fino al Passo 5.

---

## Passo 2. Aprire un terminale

Il terminale è una finestra in cui si scrivono comandi.

* **Windows:** premi il tasto **Start**, scrivi **PowerShell** e apri **Windows PowerShell**
  (oppure **Terminale** su Windows 11).
* **macOS:** premi **Cmd + Spazio**, scrivi **Terminale** e premi Invio.
* **Linux:** premi **Ctrl + Alt + T**.

Controlla che Python funzioni:

```bash
python --version
```

Dovresti vedere qualcosa come `Python 3.12.7`. Se compare un errore, vai a
[Risoluzione dei problemi](#risoluzione-dei-problemi).

---

## Passo 3. Scaricare il programma

Scegli **una** delle due possibilità.

### Opzione A: scaricare uno ZIP (la più semplice)

1. Apri <https://github.com/badoglio/Adaptive-PPG-Time-Windows>.
2. Fai clic sul pulsante verde **"Code"**, poi su **"Download ZIP"**.
3. Estrai il file ZIP in una cartella facile da ritrovare, per esempio `Documenti`.
   Su Windows: clic destro sul file → **"Estrai tutto…"**.
4. Ora hai una cartella chiamata **`Adaptive-PPG-Time-Windows-main`**.

### Opzione B: usare Git (semplifica gli aggiornamenti)

Se Git è installato (<https://git-scm.com/downloads>), digita:

```bash
cd Documents
git clone https://github.com/badoglio/Adaptive-PPG-Time-Windows.git
```

Ora hai una cartella chiamata **`Adaptive-PPG-Time-Windows`** dentro `Documenti`.

> Nel terminale la cartella `Documenti` si chiama `Documents` anche su un sistema in italiano.

---

## Passo 4. Entrare nella cartella del programma

Nel terminale, spostati nella cartella del Passo 3 con il comando `cd` ("change directory").

* **Windows:** il modo più semplice è aprire la cartella in Esplora file, fare clic destro in un punto
  vuoto e scegliere **"Apri nel terminale"**. Altrimenti digita (adatta il percorso alla tua cartella):

  ```powershell
  cd "$HOME\Documents\Adaptive-PPG-Time-Windows-main"
  ```

* **macOS / Linux:** digita `cd ` (con uno spazio finale), trascina la cartella nella finestra del
  terminale e premi Invio. Oppure digita:

  ```bash
  cd ~/Documents/Adaptive-PPG-Time-Windows-main
  ```

Con l'Opzione B il nome della cartella non ha `-main` alla fine.

Controlla di essere nel posto giusto: il comando `dir` (Windows) o `ls` (macOS/Linux) deve elencare, tra
gli altri, `pyproject.toml`, `README.md` e le cartelle `src`, `app` e `sample_data`.

> Consiglio: se il percorso contiene spazi, mettilo tra virgolette (`"…"`).

---

## Passo 5. Creare un ambiente virtuale

Un ambiente virtuale è una cartella privata (`.venv`) in cui vengono installati il programma e le sue
librerie. Li tiene separati dal resto del computer e si può cancellare in qualsiasi momento.

```bash
python -m venv .venv
```

Richiede qualche secondo e, se va a buon fine, non stampa nulla.

---

## Passo 6. Attivare l'ambiente virtuale

Va fatto **ogni volta** che apri un nuovo terminale per usare il programma.

* **Windows (PowerShell):**

  ```powershell
  .\.venv\Scripts\Activate.ps1
  ```

  Se compare un errore che dice che *l'esecuzione di script è disabilitata nel sistema*, digita una sola
  volta questo comando, rispondi **S** (Sì) se richiesto, poi ripeti l'attivazione:

  ```powershell
  Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
  ```

* **Windows (Prompt dei comandi, `cmd`):**

  ```bat
  .venv\Scripts\activate.bat
  ```

* **macOS / Linux:**

  ```bash
  source .venv/bin/activate
  ```

Quando l'ambiente è attivo, la riga del terminale inizia con **`(.venv)`**.

---

## Passo 7. Installare il programma

Con l'ambiente attivo (vedi `(.venv)`), digita:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[app]"
```

Il secondo comando scarica e installa tutto ciò che serve al programma (numpy, scipy, pandas, plotly,
Streamlit…). Può richiedere qualche minuto. Ha finito quando compare `Successfully installed …`.

Lascia le virgolette intorno a `".[app]"`: alcuni terminali (per esempio su macOS) le richiedono.

**Componenti opzionali.** Per leggere anche file EDF/BDF e file di configurazione YAML, installa invece:

```bash
python -m pip install -e ".[app,edf,yaml]"
```

> Dopo l'installazione non spostare né rinominare la cartella del programma: l'installazione punta lì. Se
> la sposti, ripeti i Passi 5–7.

---

## Passo 8. Verificare l'installazione

```bash
adaptive-ppg --version
```

Dovresti vedere `adaptive-ppg 0.2.0` (o una versione successiva).

---

## Passo 9. Avviare la dashboard

```bash
adaptive-ppg dashboard
```

* La prima volta Streamlit può chiedere un indirizzo e-mail nel terminale. È facoltativo: premi
  semplicemente **Invio**.
* Dopo qualche secondo il browser apre la pagina **<http://localhost:8501>**. Se non si apre da solo,
  copia quell'indirizzo nel browser.
* La pagina è servita dal tuo computer: i tuoi dati non vengono caricati su internet.

Per **chiudere** la dashboard, torna al terminale e premi **Ctrl + C**. Chiudere la scheda del browser non
basta.

---

## Passo 10. La prima analisi

1. Nella barra laterale a sinistra, apri **"1. Data"** e scegli **"Sample data"** come Source.
2. In **"File"** scegli **`Sample1.CSV`** (una registrazione PPG acquisita con luce infrarossa alla
   falange distale, con il soggetto seduto).
3. Lascia le altre impostazioni come sono. L'analisi parte da sola.
4. Sfoglia le schede in alto: **Signal & beats**, **Variability & windows**, **Features**,
   **Templates**, **Benchmark** ed **Export**.
5. In **Export** puoi scaricare i risultati come ZIP di file CSV, come cartella di lavoro Excel, e la
   configurazione usata (JSON).

Per analizzare una tua registrazione, scegli **"Upload a file"** in **"1. Data"** (CSV, TXT, TSV, MAT,
EDF o BDF). Se il file non ha una colonna del tempo, togli la spunta a **"Infer the sampling rate from the
file"** e scrivi la frequenza di campionamento in Hz.

L'interfaccia della dashboard è in inglese.

---

## Passo 11 (facoltativo). Usare la riga di comando

Senza la dashboard, puoi analizzare un file direttamente dal terminale:

```bash
adaptive-ppg run sample_data/Sample1.CSV --out results/sample1
```

La cartella `results/sample1` conterrà `beats.csv`, `features_per_beat.csv`, `features_grid.csv`,
`windows.csv`, `windows_fixed_30s.csv`, `feature_dictionary.csv` e `metadata.json`. Per vedere tutte le
opzioni:

```bash
adaptive-ppg run --help
```

Il [README](../README.md) principale (in inglese) descrive il metodo e ogni opzione.

---

## Le volte successive

Non serve reinstallare nulla. Apri un terminale e:

1. entra nella cartella del programma (Passo 4);
2. attiva l'ambiente (Passo 6);
3. avvia la dashboard (Passo 9).

---

## Aggiornare a una nuova versione

* **Se hai usato Git (Opzione B):** entra nella cartella, attiva l'ambiente, poi:

  ```bash
  git pull
  python -m pip install -e ".[app]"
  ```

* **Se hai usato lo ZIP (Opzione A):** scarica il nuovo ZIP, estrailo in una nuova cartella e ripeti i
  Passi 4–8. Poi cancella la vecchia cartella, dopo aver copiato altrove i risultati che vuoi conservare.

---

## Disinstallare

Cancella la cartella del programma. L'ambiente virtuale (`.venv`) è al suo interno, quindi non resta
nient'altro. Python si può rimuovere come qualsiasi altro programma.

---

## Risoluzione dei problemi

**`python` non è riconosciuto / "command not found"**
* Windows: probabilmente Python è stato installato senza "Add python.exe to PATH". Riavvia l'installer,
  scegli **Modify** (o reinstalla) e spunta quell'opzione. Puoi anche provare `py --version`: se funziona,
  al Passo 5 scrivi `py` al posto di `python`.
* Windows: se `python` apre il **Microsoft Store**, apri *Impostazioni → App → Impostazioni app avanzate →
  Alias di esecuzione app* e disattiva le due voci "python".
* macOS / Linux: digita `python3` al posto di `python`.

**`python --version` mostra una versione precedente alla 3.10**
Installa Python 3.12 (Passo 1) e usalo per creare l'ambiente. Su Windows: `py -3.12 -m venv .venv`.

**`adaptive-ppg` non è riconosciuto**
L'ambiente virtuale non è attivo: ripeti il Passo 6 (la riga deve iniziare con `(.venv)`). In alternativa,
`python -m adaptive_ppg --version` fa la stessa cosa.

**"L'esecuzione di script è disabilitata nel sistema" (Windows)**
Vedi il Passo 6, oppure usa l'attivazione dal Prompt dei comandi (`cmd`).

**`pip install` fallisce con errori di rete o SSL**
Controlla la connessione a internet. Sulla rete di un'università o di un'azienda può servire un proxy:
chiedi al servizio informatico, poi usa `python -m pip install --proxy http://indirizzo:porta -e ".[app]"`.

**`error: externally-managed-environment` (Linux / macOS)**
L'ambiente non è attivo: attivalo (Passo 6) e ripeti il comando.

**La dashboard dice che l'indirizzo è già in uso**
C'è ancora un'altra dashboard in esecuzione. Chiudila con Ctrl + C nel suo terminale, oppure avvia questa
su un'altra porta:

```bash
adaptive-ppg dashboard --server.port 8502
```

**La pagina resta bianca o mostra un errore**
Leggi il messaggio di errore nel terminale. Controlla di aver installato con `".[app]"` (Passo 7) e di aver
lanciato il comando dalla cartella del programma.

**Ancora bloccato?**
Apri una segnalazione su <https://github.com/badoglio/Adaptive-PPG-Time-Windows/issues>. Indica il sistema
operativo, l'output di `python --version`, il comando che hai digitato e il messaggio di errore completo.
