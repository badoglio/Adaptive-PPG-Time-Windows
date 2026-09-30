# Adaptive PPG Windows

> **脈波形の非定常性に基づく可変幅アンサンブル平均による光電容積脈波（PPG）の形態解析**

[English README](README.md)

本リポジトリは、Python パッケージ（`adaptive_ppg`）と Streamlit ダッシュボードで構成されています。PPG 信号から 33 の形態特徴量を抽出します。特徴量は固定区間で平均するのではなく、**脈波形の変化の速さに応じて長さが変わるアンサンブル窓**で平均します。

* 形態が安定している間は長い窓（最大 `T_max`）でノイズを抑えます。
* 過渡期は短い窓（最小 `T_min` 拍）で時間分解能を確保します。

窓長・オーバーラップ・因果的なアンカー配置は特徴量カテゴリごとに計算します。すべての記録を共通の処理周波数にリサンプリングするため、解析は入力のサンプリング周波数に依存しません。結果は、真値が既知の合成記録と、2 つの公開データセットの実記録で検証しています。

---

## 背景

**定常性**と**時間分解能**のトレードオフは、通常、固定窓（例：30 秒、または 1996 年 HRV タスクフォースの 5 分）を選んだ時点で一度に決まってしまいます。しかし自律神経の応答は定常でも同期的でもありません。

* 迷走神経の作用は数秒以内に現れます。
* 交感神経による血管収縮は数十秒かけて進みます。
* 回復には数分かかることもあります。

単一の固定窓では、速い過渡応答がぼやけるか、長い安定期にノイズが残ります。本プロジェクトは発想を逆にし、**測定した脈波形の非定常性が窓長 `T_w` とオーバーラップ `O_w` を決める**ようにしています。

---

## 手法

```
生 PPG ─► 125 Hz へのリサンプリング ─► 極性の判定
       ─► 2 帯域フィルタ ─► 拍検出（Elgendi）＋ 接線法によるサブサンプル精度の立ち上がり点 ─► 拍ごとの SQI
       ─► 正規化拍（ピーク位置を揃えて 128 点にリサンプリング）
       ─► 形態変動性（drift_z）─► H ∈ [0, 1] ─► カテゴリごとの T_w^c, O_w
       ─► 窓ごとのアンサンブルテンプレート ─► 特徴量（＋ CI95）─► 1 Hz 特徴量グリッド
```

### 1. 前処理と信号品質（`preprocessing.py`）

* **共通の処理周波数：** 信号と補助チャンネル（ECG など）を、ポリフェーズのアンチエイリアスフィルタ（`scipy.signal.resample_poly`、有理数比のため出力周波数は厳密）で `target_fs`（既定 125 Hz）にリサンプリングします。
  * 目標から 1 % 以内の周波数はそのままです。`target_fs=None` で元の周波数を保ちます。
  * クリッピングは元の信号で検出し、そのマスクを新しい時間軸に移します。
  * アップサンプリング時も、フィルタの遮断周波数は*入力*のナイキスト周波数で制限されるため、帯域は広がりません。
  * 125 Hz を選んだ理由：30 Hz の形態帯域を余裕をもって含み、BIDMC/MIMIC と同じ周波数であり、異なる機器の記録を同じ時間軸で解析できるためです。メモリ上およびエクスポートされる信号は周波数比の分だけ小さくなります。実行時間はほとんど変わりません。処理時間の大半は固定グリッド（後述）で動くテンプレート段階が占めるためです（`Sample1.CSV`：元の 500 Hz で 0.33 秒、125 Hz で 0.32 秒）。
  * 周波数の影響：`Sample1.CSV` と `Sample2.CSV` では、125 Hz での解析と元の 500 Hz での解析の差は、特徴量ごとの中央値で拍ごとの SD の 0.10 倍です（10 秒窓では 0.03–0.05 SD）。250 Hz では 0.02 SD、100 Hz では 0.14 SD、64 Hz では 0.3 SD です。微分特徴量には 100 Hz 未満は推奨しません。
* **極性：** センサによってはフォトダイオードの生カウント（光が多いほど大きい値）を記録するため、脈波が反転します。`polarity="auto"`（既定）では一次微分の歪度の符号を使います。正立した PPG は速く立ち上がりゆっくり減衰するため、歪度は正になります。負の場合は信号を反転し、警告を記録します。`"normal"` と `"inverted"` で極性を固定できます。`raw_signal` は元の極性を保つため、DC レベルと灌流指数は影響を受けません。
* **ゼロ位相 Butterworth フィルタを 2 系統（SOS）で使います。**
  * 検出用：0.5–8 Hz。Elgendi (2013) の収縮期ピーク検出器に合わせています。
  * 形態用：0.5–30 Hz。立ち上がりと二次微分波形を保持します。
  * ナイキスト周波数を超える遮断周波数は制限され、警告が記録されます。
* **立ち上がり点：** **接線法**で求めます。最大傾き点の接線と、直前の極小を通る水平線との交点です。極小点そのものは拡張期の減衰に引きずられて動くため、接線法の方が再現性に優れます。`onset_method="minimum"` も選べます。
* **サブサンプル精度の基準点：** 接線法の立ち上がり点は小数のサンプル位置として保持し、ピークと極小は放物線補間で補正します（`onset_positions`、`peak_positions`、`end_positions`。`onsets`/`peaks`/`ends` は整数インデックスのままです）。そのため拍の時刻、IBI、時間特徴量は `1 / fs`（125 Hz で 8 ms）に量子化されません。
* **拍ごとの因果的 SQI：** 次のいずれかに当たる拍は除外します。
  * 心拍数が `[min_bpm, max_bpm]` の範囲外。
  * IBI が直近の中央値から `max_ibi_jump` を超えて変化している。
  * 逐次更新される参照拍との相関が `min_template_corr` 未満。
  * 振幅比が `amplitude_ratio_limits` の範囲外。
  * クリッピングしている。
* **拍の正規化：** 各拍を `[0, 1]` にスケーリングし、区分的に（PCHIP で）リサンプリングします。収縮期ピークは常に 128 点の 30 % の位置に来るため、変動性はタイミングではなく形状の変化を反映します。
* **検出器：** 既定はネイティブ実装の Elgendi 検出器です。`detector="biosppy"` は任意です（`pip install -e ".[biosppy]"`）。BioSPPy は `peakutils` を宣言せずに import するため、この extra で一緒にインストールします。

### 2. 形態変動性と H（`engine.py`）

既定の指標 **`drift_z`** は、拍ごとのばらつきではなく**非定常性**を測ります。

```
drift_k   = RMS( mean(拍 k-N+1..k) − mean(拍 k-2N+1..k-N) )                N = N_past（10）
drift_z_k = drift_k / ( sqrt(2/N) · ブロック内 RMS 偏差のプール値 )
```

信号が定常なら 2 つのブロック平均の差はノイズだけによるもので、**ノイズレベルにかかわらず** `drift_z ≈ 1` になります。形状が持続的に変化すると、ノイズ単位で測ったその変化量に比例して値が大きくなります。

`drift_z` は**固定スケール**で `H ∈ [0, 1]` に写像します。`drift_z = 1` で `H = 0`、`drift_z = 3` で `H = 1` です（`fixed_scale=(1, 3)`）。有効拍が `2N` 個そろうまで `H` は 0 に保たれます（ウォームアップ）。

その他の選択肢：
* 指標：`drift`、`successive`（連続する拍の RMS 差）、`dispersion`。
* データ駆動の正規化：`calibration`（最初の 60 秒の p5–p95）、`rolling`、`global`。
* 変動性のうち IBI に起因する成分の回帰除去（`ibi_correction`）。

`successive` と `calibration` の組み合わせは、この README が当初記述していた手法で、`legacy` プリセットとして残しています。これには 2 つの弱点があります。
* 主に拍ごとのノイズを測ってしまいます（合成ベンチマークでは、過渡期 0.030 に対し安静時 0.026）。
* パーセンタイル正規化の性質上、安静時の記録が `H ≈ 0.5` になります。

### 3. 窓長とオーバーラップ

各カテゴリ `c` について：

```
s(H)   = 1 / (1 + exp(k · (H − θ)))                   （シグモイド；"linear" は 1 − H）  k = 10, θ = 0.3
T_w^c  = T_min^c + (T_max^c − T_min^c) · s(H)          T_max^c [拍] = T_max_sec^c / 局所 IBI
T_w^c  = max(round(T_w^c), T_crit^c)                   有効拍数で数える
O_w    = O_max − (O_max − O_min) · s(H)                短い窓ほどオーバーラップが大きい（0.25 → 0.85）
```

* **伸長の速度制限：** 窓の伸長は `max_expand_per_beat`（1 拍あたり 1 拍）に制限され、収縮は即時です。そのため、過渡応答には即座に反応しつつ、`H` がノイジーでも振動しません。
* **配置：** 隣り合うアンカーの間隔は `round(T_w · (1 − O_w))` 拍です。
* **アンカー：** 既定は**後方（trailing）窓**で、時刻 `t` の値は `t` までの拍だけに依存する因果的な出力になります。`StreamingAdaptiveEngine` は、1 拍ずつの処理でバッチ結果を正確に再現します。`anchor="centered"` を選ぶと、オフライン用の遅延ゼロの窓になります。
* **カバー範囲：** 最後の拍は常にアンカーになるため、記録の末尾まで覆われます。

| カテゴリ | 特徴量 | `T_min` | `T_max` | `T_crit` |
| :--- | :--- | :---: | :---: | :---: |
| **macro**（17） | ピーク振幅、10/25/50/75 % 脈波幅、全体/収縮期/拡張期面積、面積比、変曲点面積比、変曲点振幅、反射指数、最大収縮期/減衰傾き、傾き比、DC レベル、灌流指数 | 3 拍 | 30 秒 | 3 |
| **time_volume**（10） | 脈波持続時間、心拍数、クレスト時間、減衰時間、デューティ比、ノッチ時間、ΔT_DVP、スティフネス指数（身長が必要）、足部/ピークまでの PAT（ECG が必要） | 3 拍 | 30 秒 | 3 |
| **derivatives**（6） | SDPTG b/a, c/a, d/a, e/a、加齢指数 (b−c−d−e)/a、(b−e)/a | 5 拍 | 60 秒 | 10 |

### 4. アンサンブルテンプレートと特徴量（`features.py`）

各窓の有効拍に対して、次の処理を行います。
1. サブサンプル精度で求めた最大傾き点で位置を揃えます。
2. 処理周波数によらず、**固定のテンプレートグリッド**（`template_fs`、既定 500 Hz）上で 3 次 B スプライン補間によりサンプリングします。
3. 足部から足部への直線ベースラインを除去します。
4. ピークを 1 に正規化して平均します。
5. 平均振幅で元のスケールに戻します。

**テンプレート SQI** は、各拍とテンプレートとの相関の平均です。テンプレート上のピーク、ノッチ、拡張期の点は放物線補間で補正します。

**微分**は Savitzky–Golay（70 ms、3 次）で、テンプレートの前後にフィルタ窓 1 つ分のマージンを付けて計算します。テンプレートグリッドが固定なので、どの入力周波数でもフィルタは同じです。固定しない場合、窓は奇数サンプルに丸められ（例：250 Hz で 76 ms）、同じ合成信号で SDPTG b/a が −1.12 から −1.59 にずれていました。

> **SDPTG の比は微分の平滑化に強く依存します。** 同じ信号で、b/a は 50 ms の窓で −0.65、70 ms で −1.12、80 ms で −1.64 です。b/a、c/a、d/a、e/a、加齢指数は、同じ `sg_window_sec` を使った解析どうしでのみ比較し、その値を報告してください。

**SDPTG の各点**は制約付きで検出します。a < b < c < d < e の順序を守り、c と d には a の 5 % 以上のプロミネンスが必要です。

**重複切痕点**は、ノッチ（極小の後に極大が続く形）が見つかればその点、見つからなければ変曲点とします。

**不確かさ：** 各窓の値には、拍ごとの値から求めた `CI95 = 1.96 · SD / sqrt(n)` が付きます。

**拍系列の特徴量：**
* DC レベル、灌流指数、PAT は窓内の拍について平均します。
* PAT は、任意の ECG チャンネルから検出した R 波を使います。
* すべての特徴量について拍ごとの値も出力されます。

**出力グリッド：** 窓の値は 1 Hz のグリッドに再サンプリングします。後方窓は 0 次ホールドでグリッドの因果性を保ち、中央窓は線形補間します。

比較用の固定窓ベースライン（例：10 秒と 30 秒、オーバーラップ 50 %）も同じコードで計算します。

---

## 検証

`synthetic.py` は真値が既知の PPG を生成します。
* **拍モデル：** 1 拍を 3 つのガウス波（収縮期波、tidal 波、拡張期反射波）で表します。
* **シナリオ：** 安静 120 秒 → 急性ストレス（τ = 4 秒）→ 緩やかな回復（τ = 25 秒）。ストレス中は次のように変化します。
  * 心拍数は 70 から 100 bpm に上がり、PAT は 250 から 190 ms に下がります。
  * 振幅は 40 % 低下します。
  * tidal 波と反射波が小さくなります。
* **ノイズ付加：**
  * 白色ノイズ（30 dB）とベースライン変動。
  * 呼吸性の振幅変調と RSA、IBI のゆらぎ。
  * 体動アーチファクトのバースト（約 1 回/分）。
  * 合成 ECG。

`benchmark.py` はパイプラインを実行し、各推定器を真値と比較して評価します。真値は、同じシナリオをノイズなしで再生成した記録の拍ごとの特徴量を、1 Hz グリッドに補間したものです。推定器は次のとおりです。
* 適応窓。
* 10/30/60 秒の固定窓。
* 拍ごとの値の 10 秒 EMA。
* 拍ごとの生の値。

評価指標：
* **NRMSE：** RMSE をストレスシナリオでの真値の SD で割った値。小さいほど良い。
* **定常 NSD：** 定常区間でのノイズ。
* **遅延：** ストレス開始後に 50 % 交差が起きるまでの遅れ。

特徴量についての中央値、**ストレスシナリオ、3 シード**（`adaptive-ppg benchmark --seeds 1 2 3`）：

| 推定器 | NRMSE | 定常 NSD | 遅延（秒） |
| :--- | :---: | :---: | :---: |
| 固定 10 秒 | **0.374** | 0.158 | **8.2** |
| **適応（既定）** | 0.414 | **0.116** | 9.3 |
| EMA 10 秒 | 0.459 | 0.261 | 11.2 |
| 拍ごと | 0.637 | 0.506 | – |
| 固定 30 秒 | 0.647 | 0.383 | 24.7 |
| 固定 60 秒 | 0.971 | 0.641 | 40.5 |

**安静シナリオ**（定常、`--scenario rest`、シード 1）の NRMSE：固定 60 秒 0.079、EMA 10 秒 0.081、固定 30 秒 0.091、**適応 0.093**、固定 10 秒 0.116。

**解釈：** 両方の状況で良好な固定窓はありません。
* 10 秒窓は過渡応答に追従できますが、安静時には最もノイジーです。
* 60 秒窓は安静時に最良、ストレス時に最悪です。

適応窓は、それぞれの状況で最良の固定窓に近い性能を示し、窓長を事前に選ぶ必要がありません。ストレスプロトコル中の定常ノイズは最小です。ただし、どちらの状況でも勝者では**なく**、2 つのシナリオの平均 NRMSE（0.254）は固定 10 秒窓（0.245）に近いものの、それを下回りません。利点は定常区間にあります。10 秒窓よりノイズが小さく、その代わり遅延が約 1 秒大きくなります。

**壊れやすい特徴量：** 短い（3 拍の）窓では、壊れやすい基準点のノイズが増えます。最大の例はサブサンプル精度の基準点で解消しました。ΔT_DVP の NRMSE は、現在は適応窓で 0.46、固定 10 秒で 0.59 です（整数の基準点では 0.86 と 0.39 でした）。変曲点面積比はどの推定器でも NRMSE > 1 であり、解釈には注意が必要です。壊れやすい特徴量では、そのカテゴリの `T_crit` を大きくしてください。

### サンプリング周波数への非依存性

`tests/test_sampling.py` では特に、同じ 500 Hz の記録を元の周波数で解析しても 125 Hz で解析しても、固定窓の特徴量が一致することを確認しています。合成の安静シナリオでは、125 Hz での定常 NSD は、整数の基準点での 0.104 から、サブサンプル精度の基準点と固定テンプレートグリッドにより 0.036 に下がりました（元の 500 Hz と 1000 Hz での解析ではそれぞれ 0.041 と 0.039）。拍ごとの b/a の中央値は、1000 Hz から 100 Hz までのどの処理周波数でも −1.12 から −1.17 の間に収まります。

### 実記録

Pulse Transit Time PPG Dataset の 2 つの抜粋（`Sample1.CSV`、`Sample2.CSV`、500 Hz、末節骨での赤外光、座位。[サンプルデータ](#サンプルデータ)参照）と BIDMC のレコード 01 を既定の設定で解析しました。`tests/test_real_data.py` は、結果が妥当な範囲に収まることを確認します。

| | `Sample1.CSV` | `Sample2.CSV` | `bidmc_01` |
| :--- | :---: | :---: | :---: |
| 周波数（入力 → 処理） | 500 → 125 Hz | 500 → 125 Hz | 125 Hz |
| 極性（微分の歪度） | 反転（−2.45） | 反転（−2.39） | 正立 |
| 拍数（有効） | 335（100 %） | 277（99.6 %） | 718（98.2 %） |
| 心拍数（中央値） | 71 bpm | 74 bpm | 91 bpm |
| drift_z の中央値 / H > θ の拍 | 1.18 / 26 % | 1.06 / 21 % | 0.83 / 7 % |
| macro 窓、中央値 [p10–p90] | 17 [4–33] 拍 | 18 [5–34] 拍 | 35 拍 |
| derivatives 窓、中央値 | 24 拍 | 26 拍 | 59 拍 |
| テンプレート SQI、macro / derivatives | 0.985 / 0.958 | 0.989 / 0.982 | 0.993 |
| SDPTG b/a | −0.72 | −0.70 | −0.87 |

* PTT の 2 つの抜粋はどちらも反転して記録されており、自動の極性判定がこれを検出しました。
* drift_z は 3 つの記録すべてで 1 に近く、形態が安定しているときの期待どおりです。PTT の抜粋は BIDMC より H > θ の拍が多いため、窓が短くなります。どちらも座位の安静時に記録されたので、体動ではこの差を説明できません。センサ（指先での MAX30101 の生カウントと、臨床モニタの処理済み pleth チャンネルの違い）や接触圧の方が原因として考えやすいですが、検証はしていません。
* 1 Hz グリッド上では、SDPTG 特徴量の粗さは適応窓の方が固定 10 秒窓の約半分で（Sample1 の b/a：1 秒あたり拍ごとの SD の 0.09 倍に対し 0.21 倍）、CI95 は同程度かそれ以下です。ノッチ時間と反射指数は、短い macro/time_volume 窓に従うため、Sample1 では適応窓の方が粗くなります。

感度解析には `adaptive-ppg sweep` を使います。既定値（`drift_z`、`fixed_scale=(1, 3)`、`θ = 0.3`、`N_past = 10`、`T_min = 3`）は、ストレス 3 シードと安静 1 シードでのスイープから決めました。

---

## インストール

Python 3.10 以上を推奨します（コアパッケージは 3.9 にも対応）。Python に不慣れな方向けの手順書：[日本語](docs/INSTALL_ja.md) · [English](docs/INSTALL_en.md) · [Italiano](docs/INSTALL_it.md)。

```bash
git clone https://github.com/badoglio/Adaptive-PPG-Time-Windows.git
cd Adaptive-PPG-Time-Windows
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[app,edf,yaml]"   # pytest/ruff には "dev"、別の検出器には "biosppy" を追加
```

コアパッケージの依存は numpy、scipy、pandas、plotly だけです。extras で次のものが追加されます。
* `edf`：pyedflib（EDF/BDF ファイル用）。
* `app`：streamlit と openpyxl。
* `yaml`：YAML 形式の設定ファイル。
* `biosppy`：BioSPPy 検出器（peakutils も一緒にインストール）。
* `dev`：pytest と ruff。

## 使い方

### ダッシュボード

```bash
adaptive-ppg dashboard             # または：streamlit run app/streamlit_app.py
```

`adaptive-ppg dashboard` はソースのチェックアウト（editable インストール）で動作し、`app/streamlit_app.py` を探します。後ろに続く引数は Streamlit に渡されます。

**データソース：**
* 同梱の `sample_data`。
* 手元の CSV/TSV/TXT/MAT/EDF/BDF ファイル。サンプリング周波数は時間列やヘッダから推定するか、手入力します。CSV の区切り文字（`,` `;` タブ）、小数点のカンマ、UTF-8 BOM は自動で検出します。列を空欄にすると、時間列以外の最初の列を使います。
* 真値付きの合成デモ。

**サイドバー：**
* プリセット。
* 変動性指標、正規化、θ、k、N_past。
* 窓のアンカー。
* カテゴリごとの `T_min`/`T_max`/`T_crit`。
* フィルタ帯域、立ち上がり点の手法、検出器。
* 共通の処理周波数へのリサンプリングと、信号の極性（auto/normal/inverted）。
* 身長。
* 固定窓ベースライン。

**タブ：**
* 信号と拍（除外拍とアーチファクトを網掛け表示）。
* 変動性、H、`T_w`、`O_w`。
* 特徴量（CI95 付きの適応窓を固定窓・拍ごとの値と比較。デモでは真値も表示）。
* アンサンブルテンプレート。
* ベンチマーク。
* エクスポート（CSV と `metadata.json` の ZIP、Excel ブック、設定 JSON）。

### コマンドライン

```bash
adaptive-ppg run sample_data/bidmc_01_Signals.csv --signal-column PLETH --ecg-column II --out results/
adaptive-ppg run "data/*.edf" --channel Pleth --preset acute_stress --fixed 10 30 --templates --out results/
adaptive-ppg run my_recording.csv --signal-column PPG --fs 500 --out results/mine   # 時間列がない場合は fs を指定
adaptive-ppg run sample_data/Sample1.CSV --target-fs 250 --out results/sample1   # --target-fs 0：元の周波数
adaptive-ppg init-config my_config.json --preset rest        # 編集してから --config で指定
adaptive-ppg benchmark --seeds 1 2 3 --out bench/
adaptive-ppg sweep --param engine.theta 0.2 0.3 0.5 --param engine.N_past 6 10 --with-rest --out sweep.csv
```

バッチ実行では、記録ごとのフォルダと `batch_summary.csv` が出力されます。すべてのエクスポートの `metadata.json` には次の情報が含まれます。
* 全設定。
* パッケージのバージョンと git リビジョン。
* 入力のサンプリング周波数とその取得元、および処理周波数。
* 検出した極性と、その根拠となった微分の歪度。
* 警告。

### Python

```python
from adaptive_ppg.io import BioSignalLoader
from adaptive_ppg.config import preset
from adaptive_ppg.pipeline import run

sd = BioSignalLoader.auto_load("sample_data/bidmc_01_Signals.csv", signal_column="PLETH", ecg_column="II")
res = run(sd, preset("default"), fixed_window_secs=[30])
res.features.grid            # 1 Hz の特徴量グリッド（適応窓）
res.window_table()           # 全適応窓とその特徴量・CI95
res.save("results/bidmc01", include_templates=True)
```

リアルタイム用途では、`adaptive_ppg.engine.StreamingAdaptiveEngine.push(...)` に 1 拍ずつ渡すと、閉じた窓が順に返されます。

### 設定

すべてのパラメータは、シリアライズ可能な 1 つの `PipelineConfig` にまとまっています。セクションは `preprocessing`、`engine`（`categories` を含む）、`features` です。部分的な JSON/YAML ファイルは既定値にマージされ、未知のキーはエラーになり、`null` を指定したカテゴリは無効になります。

サンプリング周波数に関係するパラメータ：

| パラメータ | 既定値 | 意味 |
| :--- | :---: | :--- |
| `preprocessing.target_fs` | 125.0 | 処理周波数（Hz、20 以上）。`null` で元の周波数を保つ |
| `preprocessing.polarity` | `"auto"` | `"auto"`、`"normal"`、`"inverted"` |
| `features.template_fs` | 500.0 | テンプレートグリッド（Hz、50 以上）。`null` で処理周波数を使う |
| `features.sg_window_sec` | 0.07 | 微分の Savitzky–Golay 窓（SDPTG の結果とともに報告すること） |

プリセット：
* `default`。
* `rest`：長めの窓、θ = 0.5。
* `acute_stress`：短めの窓、θ = 0.2。
* `wearable_64hz`：狭い形態帯域（15 Hz）、長い微分平滑化（110 ms）、derivatives の大きい `T_crit`。
* `legacy`：当初の successive/calibration 方式。

---

## パッケージ構成

```
src/adaptive_ppg/
  io.py             CSV/TSV/TXT（書式の自動判定）、MAT（v5–v7.2）、EDF/BDF の読み込み；頑健な fs 推定
  preprocessing.py  リサンプリング、極性、フィルタ、Elgendi 検出、サブサンプル精度の基準点、SQI、
                    拍の正規化
  engine.py         変動性指標、H、カテゴリごとのスケジュール；バッチ/ストリーミングエンジン
  features.py       アンサンブルテンプレート、33 特徴量、PAT、固定窓ベースライン、出力グリッド
  synthetic.py      真値付き合成 PPG/ECG
  benchmark.py      真値、評価指標、ベンチマーク、パラメータスイープ
  pipeline.py       run() と PipelineResult（表、ZIP/Excel/CSV 出力、来歴情報）
  visualization.py  ダッシュボード用 Plotly 図
  config.py         PipelineConfig、検証、プリセット
  cli.py            コマンドラインインターフェース
  _typing.py        配列の型エイリアス（FloatArray、IntArray、BoolArray）。パッケージは py.typed を含む
app/streamlit_app.py  Streamlit ダッシュボード
sample_data/        サンプル記録とそのライセンス（sample_data/README.md）
tests/              pytest テスト（実行：pytest -q）
```

---

## サンプルデータ

記録は第三者のデータで、それぞれ独自のライセンスに従います。[`sample_data/README.md`](sample_data/README.md) を参照してください。

* **`sample_data/Sample1.CSV`、`Sample2.CSV`**：**Pulse Transit Time PPG Dataset**（PhysioNet、v1.1.0）の抜粋（500 Hz で 282 秒と 226 秒）。Open Database License（ODbL 1.0）で公開されています。列は `Time (s)` と `Pleth` で、区切り文字は `;`、小数点はカンマです。
  * 計測条件：赤外（IR）光、末節骨（指先）、座位の安静時。
  * `Pleth` は MAX30101 の生カウントなので脈波が反転しています。自動の極性判定がこれを反転します。
  * サンプルの約 23 % が直前の値とまったく同じです。重複は均等に分布し、ほとんどが 1 サンプルだけなので、500 Hz の出力クロックへの再同期によるものと考えられます。リサンプリング後に目に見える影響はありません。
  * 引用：Mehrgardt, P., Khushi, M., Poon, S., & Withana, A. (2022). *Pulse Transit Time PPG Dataset* (version 1.1.0). PhysioNet. https://doi.org/10.13026/jpan-6n92 および Goldberger, A. L. et al., "PhysioBank, PhysioToolkit, and PhysioNet", *Circulation* 101(23):e215–e220, 2000.
* **`sample_data/bidmc_01_Signals.csv`**：PhysioNet の **BIDMC PPG and Respiration Dataset** のレコード 01（8 分、125 Hz）。列は `Time [s]`、`RESP`、`PLETH`（PPG）、`V`、`AVR`、`II`（ECG）です。
  * 時間列は小数 2–3 桁に丸められています。そのためローダは、時刻をサンプル番号に対して最小二乗フィットして fs を推定します（差分の中央値だと 100 Hz になってしまいます）。
  * このレコードの PAT は約 0.53 秒です。生理的な PAT より長く、MIMIC モニタの pleth チャンネルに知られている信号処理遅延と一致します。BIDMC/MIMIC から得た PAT の絶対値は慎重に扱ってください。
  * 引用：Pimentel, M. A. F., Charlton, P. H., & Clifton, D. A. (2016). *BIDMC PPG and Respiration Dataset* (v1.0.0). PhysioNet. https://doi.org/10.13026/C2GG69 および Pimentel et al., "Toward a robust estimation of respiratory rate from pulse oximeters", *IEEE Trans. Biomed. Eng.* 64(8):1914–1923, 2017.

## 開発

```bash
pip install -e ".[dev,edf,yaml,app]"
ruff check src tests app
pytest -q
```

GitHub Actions が Linux と Windows で ruff とテストを実行します（`.github/workflows/ci.yml`）。バージョンごとの変更点は [`CHANGELOG.md`](CHANGELOG.md)（英語）にまとめています。

## ライセンス

コードは **GNU General Public License v3.0 以降**（`GPL-3.0-or-later`、[`LICENSE`](LICENSE) 参照）で公開しています。サンプル記録はそれぞれ独自のライセンス（BIDMC は ODC-By 1.0、Pulse Transit Time PPG Dataset は ODbL 1.0）に従い、[`sample_data/README.md`](sample_data/README.md) に一覧があります。
