# インストールガイド（手順ごと）

[English](INSTALL_en.md) · [Italiano](INSTALL_it.md) · [日本語](INSTALL_ja.md)

このガイドは、Python のプログラムをインストールしたことがない方向けです。所要時間は約 15 分です。
インターネット接続と、約 1 GB の空きディスク容量が必要です。

最後には次のものが使えるようになります。

* **ダッシュボード**：自分のコンピュータ上で動き、PPG 記録を解析する Web ページ。
* **`adaptive-ppg` コマンド**：ターミナルからファイルを解析するためのコマンド。

入力するコマンドは灰色の枠内に示します。ターミナルに入力（またはコピー＆ペースト）して **Enter** を押して
ください。OS によって手順が異なる箇所は、お使いの OS の部分だけを実行してください。

---

## ステップ 1. Python をインストールする

このプログラムには **Python 3.10 以降** が必要です。動作確認に使っている **Python 3.12** をおすすめします。

### Windows

1. <https://www.python.org/downloads/windows/> を開き、Python 3.12 の **Windows installer (64-bit)** を
   ダウンロードします。
2. ダウンロードしたファイルを実行します。
3. **重要：** 最初の画面の下にある **「Add python.exe to PATH」** にチェックを入れてから、
   **「Install Now」** をクリックします。
4. 「Setup was successful」と表示されたら **Close** をクリックします。

### macOS

1. <https://www.python.org/downloads/macos/> を開き、Python 3.12 の
   **macOS 64-bit universal2 installer** をダウンロードします。
2. `.pkg` ファイルを開き、インストーラの指示に従います。
3. macOS ではコマンド名が `python` ではなく **`python3`** です。このガイドで `python` と書かれている箇所は
   `python3` と入力してください（ステップ 5 まで。仮想環境の中では `python` も使えます）。

### Linux（Ubuntu / Debian）

通常、Python はすでにインストールされています。仮想環境用のモジュールもインストールします。

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip
```

macOS と同様に、ステップ 5 までは `python` の代わりに `python3` と入力してください。

---

## ステップ 2. ターミナルを開く

ターミナルは、コマンドを入力するためのウィンドウです。

* **Windows：** **スタート** キーを押し、**PowerShell** と入力して **Windows PowerShell** を開きます
  （Windows 11 では **ターミナル** でも可）。
* **macOS：** **Cmd + Space** を押し、**ターミナル**（Terminal）と入力して Enter を押します。
* **Linux：** **Ctrl + Alt + T** を押します。

Python が動くか確認します。

```bash
python --version
```

`Python 3.12.7` のように表示されれば成功です。エラーが出た場合は [トラブルシューティング](#トラブルシューティング)
を参照してください。

---

## ステップ 3. プログラムをダウンロードする

次の 2 つの方法の **どちらか一方** を選びます。

### 方法 A：ZIP をダウンロードする（いちばん簡単）

1. <https://github.com/badoglio/Adaptive-PPG-Time-Windows> を開きます。
2. 緑色の **「Code」** ボタンをクリックし、**「Download ZIP」** をクリックします。
3. ZIP ファイルを、見つけやすいフォルダ（例：`ドキュメント`）に展開します。
   Windows の場合：ファイルを右クリック →**「すべて展開…」**。
4. **`Adaptive-PPG-Time-Windows-main`** という名前のフォルダができます。

### 方法 B：Git を使う（更新が簡単になります）

Git がインストールされている場合（<https://git-scm.com/downloads>）、次のように入力します。

```bash
cd Documents
git clone https://github.com/badoglio/Adaptive-PPG-Time-Windows.git
```

`Documents` の中に **`Adaptive-PPG-Time-Windows`** というフォルダができます。

> ターミナル上では、「ドキュメント」フォルダは日本語環境でも `Documents` という名前です。

---

## ステップ 4. プログラムのフォルダに移動する

ターミナルで、`cd` コマンド（change directory）を使ってステップ 3 のフォルダに移動します。

* **Windows：** いちばん簡単なのは、エクスプローラーでフォルダを開き、何もない所を右クリックして
  **「ターミナルで開く」** を選ぶ方法です。または次のように入力します（パスは自分のフォルダに合わせて
  変更してください）。

  ```powershell
  cd "$HOME\Documents\Adaptive-PPG-Time-Windows-main"
  ```

* **macOS / Linux：** `cd `（最後にスペース）と入力し、フォルダをターミナルのウィンドウにドラッグして
  Enter を押します。または次のように入力します。

  ```bash
  cd ~/Documents/Adaptive-PPG-Time-Windows-main
  ```

方法 B の場合、フォルダ名の末尾に `-main` は付きません。

正しい場所にいるか確認します。`dir`（Windows）または `ls`（macOS/Linux）を実行すると、
`pyproject.toml`、`README.md`、フォルダ `src`、`app`、`sample_data` などが表示されるはずです。

> ヒント：パスにスペースが含まれる場合は、引用符（`"…"`）で囲んでください。

---

## ステップ 5. 仮想環境を作る

仮想環境とは、プログラムとそのライブラリをインストールする専用のフォルダ（`.venv`）です。
コンピュータの他の部分と分けて管理でき、いつでも削除できます。

```bash
python -m venv .venv
```

数秒で終わり、成功した場合は何も表示されません。

---

## ステップ 6. 仮想環境を有効にする

プログラムを使うために新しいターミナルを開くたびに、**毎回** 行います。

* **Windows（PowerShell）：**

  ```powershell
  .\.venv\Scripts\Activate.ps1
  ```

  「このシステムではスクリプトの実行が無効になっている」というエラーが出た場合は、次のコマンドを一度だけ
  入力し、確認されたら **Y**（はい）と答えてから、もう一度有効化してください。

  ```powershell
  Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
  ```

* **Windows（コマンド プロンプト、`cmd`）：**

  ```bat
  .venv\Scripts\activate.bat
  ```

* **macOS / Linux：**

  ```bash
  source .venv/bin/activate
  ```

仮想環境が有効になると、ターミナルの行頭に **`(.venv)`** と表示されます。

---

## ステップ 7. プログラムをインストールする

仮想環境が有効な状態（`(.venv)` が表示されている状態）で、次のように入力します。

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[app]"
```

2 つ目のコマンドで、プログラムに必要なもの（numpy、scipy、pandas、plotly、Streamlit など）がすべて
ダウンロード・インストールされます。数分かかることがあります。`Successfully installed …` と表示されたら
完了です。

`".[app]"` の引用符は省略しないでください。一部のターミナル（macOS など）では必要です。

**オプション機能。** EDF/BDF ファイルや YAML 形式の設定ファイルも読み込みたい場合は、代わりに次を
インストールします。

```bash
python -m pip install -e ".[app,edf,yaml]"
```

> インストール後は、プログラムのフォルダを移動したり名前を変えたりしないでください。インストールは
> そのフォルダを参照しています。移動した場合は、ステップ 5〜7 をやり直してください。

---

## ステップ 8. インストールを確認する

```bash
adaptive-ppg --version
```

`adaptive-ppg 0.2.0`（またはそれ以降のバージョン）と表示されれば成功です。

---

## ステップ 9. ダッシュボードを起動する

```bash
adaptive-ppg dashboard
```

* 初回は、Streamlit がターミナルでメールアドレスを尋ねることがあります。入力は任意なので、そのまま
  **Enter** を押してください。
* 数秒後、ブラウザで **<http://localhost:8501>** が開きます。自動で開かない場合は、このアドレスを
  ブラウザにコピーしてください。
* このページは自分のコンピュータから配信されています。データがインターネットにアップロードされることは
  ありません。

ダッシュボードを **終了** するには、ターミナルに戻って **Ctrl + C** を押します。ブラウザのタブを閉じる
だけでは終了しません。

---

## ステップ 10. はじめての解析

1. 左のサイドバーで **「1. Data」** を開き、Source に **「Sample data」** を選びます。
2. **「File」** で **`Sample1.CSV`** を選びます（座位の被験者の末節骨（指先）で、赤外光により計測した
   PPG 記録です）。
3. ほかの設定はそのままにします。解析は自動で始まります。
4. 上部のタブを順に見てください：**Signal & beats**、**Variability & windows**、**Features**、
   **Templates**、**Benchmark**、**Export**。
5. **Export** では、結果を CSV ファイルの ZIP や Excel ブックとして、また使用した設定（JSON）を
   ダウンロードできます。

自分の記録を解析するには、**「1. Data」** で **「Upload a file」** を選びます（CSV、TXT、TSV、MAT、EDF、
BDF）。ファイルに時間列がない場合は、**「Infer the sampling rate from the file」** のチェックを外し、
サンプリング周波数（Hz）を入力してください。

ダッシュボードの画面は英語です。

---

## ステップ 11（任意）. コマンドラインで使う

ダッシュボードを使わずに、ターミナルから直接ファイルを解析できます。

```bash
adaptive-ppg run sample_data/Sample1.CSV --out results/sample1
```

`results/sample1` フォルダに `beats.csv`、`features_per_beat.csv`、`features_grid.csv`、`windows.csv`、
`windows_fixed_30s.csv`、`feature_dictionary.csv`、`metadata.json` が作られます。すべてのオプションを
見るには次を実行します。

```bash
adaptive-ppg run --help
```

手法とすべてのオプションは [README（日本語）](../README_JA.md) で説明しています。

---

## 次回からの使い方

再インストールは不要です。ターミナルを開いて、次の順に実行します。

1. プログラムのフォルダに移動する（ステップ 4）。
2. 仮想環境を有効にする（ステップ 6）。
3. ダッシュボードを起動する（ステップ 9）。

---

## 新しいバージョンに更新する

* **Git を使った場合（方法 B）：** フォルダに移動し、仮想環境を有効にしてから次を実行します。

  ```bash
  git pull
  python -m pip install -e ".[app]"
  ```

* **ZIP を使った場合（方法 A）：** 新しい ZIP をダウンロードして新しいフォルダに展開し、ステップ 4〜8 を
  繰り返します。残したい結果を別の場所にコピーしてから、古いフォルダを削除してください。

---

## アンインストール

プログラムのフォルダを削除します。仮想環境（`.venv`）はその中にあるので、ほかには何も残りません。
Python 本体は、ほかのプログラムと同じ方法で削除できます。

---

## トラブルシューティング

**`python` が認識されない／「command not found」**
* Windows：「Add python.exe to PATH」にチェックを入れずにインストールした可能性があります。インストーラを
  もう一度実行し、**Modify** を選ぶ（または再インストールする）と、そのオプションにチェックを入れられます。
  `py --version` も試してください。これが動く場合は、ステップ 5 で `python` の代わりに `py` と入力します。
* Windows：`python` と入力すると **Microsoft Store** が開く場合は、*設定 → アプリ → アプリの詳細設定 →
  アプリ実行エイリアス* を開き、「python」の 2 つの項目をオフにします。
* macOS / Linux：`python` の代わりに `python3` と入力してください。

**`python --version` の結果が 3.10 より古い**
Python 3.12 をインストールし（ステップ 1）、それを使って仮想環境を作ります。Windows の場合：
`py -3.12 -m venv .venv`。

**`adaptive-ppg` が認識されない**
仮想環境が有効になっていません。ステップ 6 をやり直してください（行頭に `(.venv)` が表示されるはずです）。
代わりに `python -m adaptive_ppg --version` でも同じことができます。

**「このシステムではスクリプトの実行が無効になっている」（Windows）**
ステップ 6 を参照するか、コマンド プロンプト（`cmd`）での有効化を使ってください。

**`pip install` がネットワークエラーや SSL エラーで失敗する**
インターネット接続を確認してください。大学や企業のネットワークではプロキシが必要な場合があります。
情報システム部門に確認のうえ、`python -m pip install --proxy http://アドレス:ポート -e ".[app]"` を
使ってください。

**`error: externally-managed-environment`（Linux / macOS）**
仮想環境が有効になっていません。有効にしてから（ステップ 6）、コマンドをやり直してください。

**ダッシュボードで「アドレスがすでに使用されている」と表示される**
別のダッシュボードがまだ動いています。そのターミナルで Ctrl + C を押して終了するか、別のポートで
起動してください。

```bash
adaptive-ppg dashboard --server.port 8502
```

**ページが真っ白のまま、またはエラーが表示される**
ターミナルに表示されたエラーメッセージを確認してください。`".[app]"` を付けてインストールしたか
（ステップ 7）、プログラムのフォルダでコマンドを実行したかを確認してください。

**それでも解決しない場合**
<https://github.com/badoglio/Adaptive-PPG-Time-Windows/issues> で issue を作成してください。OS、
`python --version` の出力、入力したコマンド、エラーメッセージの全文を添えてください。
