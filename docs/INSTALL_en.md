# Installation guide, step by step

[English](INSTALL_en.md) · [Italiano](INSTALL_it.md) · [日本語](INSTALL_ja.md)

This guide is for people who have never installed a Python program. It takes about 15 minutes. You need an
internet connection and about 1 GB of free disk space.

At the end you will have:

* the **dashboard**, a web page that runs on your computer and analyses PPG recordings;
* the **`adaptive-ppg` command**, for analysing files from the terminal.

Commands to type are shown in grey boxes. Type them (or copy and paste them) into the terminal and press
**Enter**. Wherever the guide differs by operating system, follow only the part for yours.

---

## Step 1. Install Python

The program needs **Python 3.10 or later**. We recommend **Python 3.12**, which is the version we test.

### Windows

1. Open <https://www.python.org/downloads/windows/> and download the **Windows installer (64-bit)** for
   Python 3.12.
2. Run the downloaded file.
3. **Important:** on the first screen, tick **"Add python.exe to PATH"** at the bottom, then click
   **"Install Now"**.
4. When the installer says "Setup was successful", click **Close**.

### macOS

1. Open <https://www.python.org/downloads/macos/> and download the **macOS 64-bit universal2 installer**
   for Python 3.12.
2. Open the `.pkg` file and follow the installer.
3. On macOS the command is called **`python3`**, not `python`. Wherever this guide says `python`, type
   `python3` instead (only until Step 5: inside the virtual environment `python` also works).

### Linux (Ubuntu / Debian)

Python is usually already installed. Also install the virtual-environment module:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip
```

As on macOS, type `python3` instead of `python` until Step 5.

---

## Step 2. Open a terminal

The terminal is a window where you type commands.

* **Windows:** press the **Start** key, type **PowerShell** and open **Windows PowerShell**
  (or **Terminal** on Windows 11).
* **macOS:** press **Cmd + Space**, type **Terminal** and press Enter.
* **Linux:** press **Ctrl + Alt + T**.

Check that Python works:

```bash
python --version
```

You should see something like `Python 3.12.7`. If you see an error, go to
[Troubleshooting](#troubleshooting).

---

## Step 3. Download the program

Choose **one** of the two options.

### Option A: download a ZIP (simplest)

1. Open <https://github.com/badoglio/Adaptive-PPG-Time-Windows>.
2. Click the green **"Code"** button, then **"Download ZIP"**.
3. Extract the ZIP file into a folder you can find easily, for example `Documents`.
   On Windows: right-click the file → **"Extract All…"**.
4. You now have a folder called **`Adaptive-PPG-Time-Windows-main`**.

### Option B: use Git (makes updating easier)

If Git is installed (<https://git-scm.com/downloads>), type:

```bash
cd Documents
git clone https://github.com/badoglio/Adaptive-PPG-Time-Windows.git
```

You now have a folder called **`Adaptive-PPG-Time-Windows`** inside `Documents`.

---

## Step 4. Go into the program folder

In the terminal, move into the folder from Step 3 with the `cd` command ("change directory").

* **Windows:** the easiest way is to open the folder in File Explorer, right-click an empty area and
  choose **"Open in Terminal"**. Otherwise type (adapt the path to your folder):

  ```powershell
  cd "$HOME\Documents\Adaptive-PPG-Time-Windows-main"
  ```

* **macOS / Linux:** type `cd ` (with a space at the end), drag the folder into the terminal window and
  press Enter. Or type:

  ```bash
  cd ~/Documents/Adaptive-PPG-Time-Windows-main
  ```

With Option B the folder name has no `-main` at the end.

Check that you are in the right place: the command `dir` (Windows) or `ls` (macOS/Linux) must list, among
others, `pyproject.toml`, `README.md` and the folders `src`, `app` and `sample_data`.

> Tip: put the path in quotes (`"…"`) if it contains spaces.

---

## Step 5. Create a virtual environment

A virtual environment is a private folder (`.venv`) where the program and its libraries are installed.
It keeps them separate from the rest of your computer, and you can delete it at any time.

```bash
python -m venv .venv
```

This takes a few seconds and prints nothing when it succeeds.

---

## Step 6. Activate the virtual environment

You have to do this **every time** you open a new terminal to use the program.

* **Windows (PowerShell):**

  ```powershell
  .\.venv\Scripts\Activate.ps1
  ```

  If you see an error saying that *running scripts is disabled on this system*, type this command once and
  answer **Y** (Yes) if asked, then repeat the activation:

  ```powershell
  Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
  ```

* **Windows (Command Prompt, `cmd`):**

  ```bat
  .venv\Scripts\activate.bat
  ```

* **macOS / Linux:**

  ```bash
  source .venv/bin/activate
  ```

When the environment is active, the line in the terminal starts with **`(.venv)`**.

---

## Step 7. Install the program

With the environment active (you can see `(.venv)`), type:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[app]"
```

The second command downloads and installs everything the program needs (numpy, scipy, pandas, plotly,
Streamlit…). It can take a few minutes. It is finished when you see `Successfully installed …`.

Keep the quotes around `".[app]"`: some terminals (for example macOS) need them.

**Optional extras.** To also read EDF/BDF files and YAML configuration files, install instead:

```bash
python -m pip install -e ".[app,edf,yaml]"
```

> Do not move or rename the program folder after installing: the installation points to it. If you move
> it, repeat Steps 5–7.

---

## Step 8. Check the installation

```bash
adaptive-ppg --version
```

You should see `adaptive-ppg 0.2.0` (or a later version).

---

## Step 9. Start the dashboard

```bash
adaptive-ppg dashboard
```

* The first time, Streamlit may ask for an e-mail address in the terminal. It is optional: just press
  **Enter**.
* After a few seconds your browser opens the page **<http://localhost:8501>**. If it does not open by
  itself, copy that address into the browser.
* The page is served by your own computer: your data are not uploaded to the internet.

To **stop** the dashboard, go back to the terminal and press **Ctrl + C**. Closing the browser tab is
not enough.

---

## Step 10. Your first analysis

1. In the left sidebar, open **"1. Data"** and choose **"Sample data"** as the Source.
2. In **"File"**, choose **`Sample1.CSV`** (a PPG recording taken with infrared light at the fingertip, in
   a seated subject).
3. Leave the other settings as they are. The analysis starts on its own.
4. Browse the tabs at the top: **Signal & beats**, **Variability & windows**, **Features**,
   **Templates**, **Benchmark** and **Export**.
5. In **Export** you can download the results as a ZIP of CSV files, as an Excel workbook, and the
   configuration you used (JSON).

To analyse your own recording, choose **"Upload a file"** in **"1. Data"** (CSV, TXT, TSV, MAT, EDF or
BDF). If the file has no time column, untick **"Infer the sampling rate from the file"** and type the
sampling rate in Hz.

---

## Step 11 (optional). Use the command line

Without the dashboard, you can analyse a file directly from the terminal:

```bash
adaptive-ppg run sample_data/Sample1.CSV --out results/sample1
```

The `results/sample1` folder will contain `beats.csv`, `features_per_beat.csv`, `features_grid.csv`,
`windows.csv`, `windows_fixed_30s.csv`, `feature_dictionary.csv` and `metadata.json`. To see all the
options:

```bash
adaptive-ppg run --help
```

The main [README](../README.md) describes the method and every option.

---

## Next time

You do not need to install anything again. Open a terminal and:

1. go into the program folder (Step 4);
2. activate the environment (Step 6);
3. start the dashboard (Step 9).

---

## Updating to a new version

* **If you used Git (Option B):** go into the folder, activate the environment, then:

  ```bash
  git pull
  python -m pip install -e ".[app]"
  ```

* **If you used the ZIP (Option A):** download the new ZIP, extract it into a new folder and repeat
  Steps 4–8. Then delete the old folder, after copying out any results you want to keep.

---

## Uninstalling

Delete the program folder. The virtual environment (`.venv`) is inside it, so nothing else is left
behind. Python itself can be removed like any other program.

---

## Troubleshooting

**`python` is not recognized / "command not found"**
* Windows: Python was probably installed without "Add python.exe to PATH". Run the installer again,
  choose **Modify** (or reinstall) and tick that option. You can also try `py --version`: if it works, type
  `py` instead of `python` in Step 5.
* Windows: if `python` opens the **Microsoft Store**, open *Settings → Apps → Advanced app settings →
  App execution aliases* and switch off the two "python" entries.
* macOS / Linux: type `python3` instead of `python`.

**`python --version` shows a version older than 3.10**
Install Python 3.12 (Step 1) and use it to create the environment. On Windows: `py -3.12 -m venv .venv`.

**`adaptive-ppg` is not recognized**
The virtual environment is not active: repeat Step 6 (the line must start with `(.venv)`). As an
alternative, `python -m adaptive_ppg --version` does the same thing.

**"Running scripts is disabled on this system" (Windows)**
See Step 6, or use the Command Prompt (`cmd`) activation instead.

**`pip install` fails with network or SSL errors**
Check your internet connection. On a university or company network a proxy may be needed: ask your IT
service, then use `python -m pip install --proxy http://address:port -e ".[app]"`.

**`error: externally-managed-environment` (Linux / macOS)**
The environment is not active: activate it (Step 6) and repeat the command.

**The dashboard says the address is already in use**
Another dashboard is still running. Close it with Ctrl + C in its terminal, or start this one on another
port:

```bash
adaptive-ppg dashboard --server.port 8502
```

**The page stays blank or shows an error**
Look at the terminal for the error message. Check that you installed with `".[app]"` (Step 7) and that you
started the command from the program folder.

**Still stuck?**
Open an issue at <https://github.com/badoglio/Adaptive-PPG-Time-Windows/issues>. Include your operating
system, the output of `python --version`, the command you typed and the full error message.
