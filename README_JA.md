# Adaptive PPG Windows 🌀

> **形態的複雑性に基づく可変時間窓と動的オーバーラップによる適応的バイオシグナル (PPG) 解析**

このリポジトリは、光電脈波（PPG）信号から形態的特徴を抽出・分析するための科学的パイプラインの完全な実装を含んでいます。従来の固定時間窓に依存する分析システムとは異なり、本フレームワークは被験者の血管状態の遷移速度にリアルタイムで自律的に適応する**可変幅アンサンブル平均化**アルゴリズムを適用します。

---

## 📖 生理学的論拠

実世界のバイオシグナル解析における主要なボトルネックは、**信号の定常性**と**時間分解能**の間の恣意的なトレードオフです。

### 1. 従来の固定窓の限界
従来の解析窓（心拍変動解析 - HRV で推奨される *1996 Task Force* の 5分間窓など）は、低周波（$0.04 \text{ Hz}$ の LF帯域）における安定したスペクトル推定を保証するために設計されました。しかし、自律神経系（ANS）の急性ストレス応答は、定常的でもサブシステム間で同調的でもありません：
*   **迷走神経（副交感神経）活動：** 1秒未満または数秒単位で瞬時に反応します（HF帯域: $0.15 - 0.40 \text{ Hz}$）。
*   **交感神経活動：** 反応の遅延（レイテンシ）が $3 - 5 \text{秒}$ あり、応答が数十秒から数分に及ぶ非常に遅いダイナミクスを示します（LF帯域: $0.04 - 0.15 \text{ Hz}$）。

これら2つの経路を単一の固定窓（30秒や5分など）の中で処理することは、急速な遷移を過度に平均化して交感神経の急激なスパイクを鈍らせ、相関のない遅い変化と速い変化を混在させることを意味します。

### 2. 技術的制約から生理学へ
人間の生理機能を硬直的なスペクトル窓に適合させるのではなく、本プロジェクトはパラダイムを反転させます。**観察される特徴のダイナミクス（瞬間的な血管状態の複雑性）がセグメント長 $T_w$ とオーバーラップ率 $O_w$ を決定します。**
*   **定常フェーズ（複雑性が低い場合）：** 窓幅は拡大します（$T_w \to T_{max}$、例: 30〜40心拍）。アンサンブル平均の対象周期が増え、ランダムノイズを強力に抑制して信号対雑音比（$SNR$）を最大化します。
*   **遷移フェーズ（複雑性が高い場合、例: 急性ストレスや血管収縮）：** 窓幅は急速に収縮します（$T_w \to T_{min}$、例: 3〜5心拍）。平均化の窓を最小限に絞り、血管状態の瞬間的な進化を遅れなく追跡します。

---

## 🧠 アルゴリズムパイプライン

「窓幅の決定に特徴量が必要であり、その特徴量の計算には窓幅が必要である」という**循環参照の罠**を回避するため、本パイプラインは**2段階適応境界フレームワーク（Two-Stage Bounded Adaptive Framework）**を採用しています：

```
[ 生PPG信号 ]
      │
      ▼ (io_loader & preprocessing)
[ セグメント化された正規化脈波 [0, 1] ]
      │
      ▼
[ 段階1: 形態的エントロピー H_morph の計算 ]  ◄── 事前非定常性インデックス
      │
      ├───────────────────────────────────────┐
      ▼ (段階2: 線形マッピング)                  ▼ (マルチスケール分離)
[ ターゲット窓 Tw とオーバーラップ Ow ]       [ 臨界感度境界 T_crit^k ]
      │                                       │
      └───────────────────┬───────────────────┘
                          ▼
            [ 実行窓幅 T_w^k ]
                          │
                          ▼ (feature_extraction)
           [ 特徴量抽出とベンチマーク ]
```

### 段階A: 瞬間複雑性の評価 ($H_{morph}[n]$)
セグメント化され、振幅 $[0,1]$ と時間幅（128サンプル）が正規化された拍動行列に対して、過去 $N_{past}$ 拍（デフォルト $10$ 拍）の**形態的エントロピー ($H_{morph}$)** を計算します。アルゴリズムは連続する拍動プロファイル間の平方平均二乗（RMS）差を評価します：
*   $H_{morph} \approx 0$ （定常状態）：連続する拍動が酷似しています。血管運動トーンは安定しています。
*   $H_{morph} \to 1$ （遷移状態）：形態が拍動ごとに急速に変化しています（血管収縮やコンプライアンスの調整を示します）。

### 段階B: 窓幅（$T_w$）とオーバーラップ（$O_w$）の写像
ターゲット窓サイズ（拍動数）とオーバーラップ率は以下のように動的に写像されます：
$$T_w[n] = T_{max} - (T_{max} - T_{min}) \cdot H_{morph}[n]$$
$$O_w[n] = O_{min} + (O_{max} - O_{min}) \cdot H_{morph}[n]$$
急速な遷移時（$H_{morph} \to 1$）には、窓が収縮し、オーバーラップが最大 $80-90\%$ まで増加して遷移をスムーズに追跡します。

### 段階C: 生理学的分離と臨界感度制限 ($T_{crit}^k$)
時間周波数の不確定性（窓が短すぎると数値微分が不安定になる現象）を解消するため、特徴量のカテゴリごとに臨界感度境界 $T_{crit}^k$ を適用します：
$$T_w^k[n] = \max\left( T_{target}[n], \; T_{crit}^k \right)$$

| 特徴量カテゴリ | 例 | 理論的 $T_{min}$ | 定常時 $T_{max}$ | 臨界境界 $T_{crit}^k$ | 生理学的意味 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **マクロ形態的** | *Pulse Slope Index*, *Pulse Width 50%*, *Area Ratio* | $3 \text{ 拍}$ | $30 \text{ 秒}$ | $\ge 3 \text{ 拍}$ | 一次的な血管運動トーンとコンプライアンス。 |
| **時間容量的** | *Systolic Peak Time*, *Pulse Transit Duration* | $3 \text{ 拍}$ | $30 \text{ 秒}$ | $\ge 3 \text{ 拍}$ | 心室駆出速度と伝播特性。 |
| **微分特徴 (VPG / SDPTG)** | 二次微分からの $b/a$ や $d/a$ 比 | $5 \text{ 拍}$ | $60 \text{ 秒}$ | $\ge 10 \text{ 拍}$ | 血管の硬化度（SDPTGは数値安定のために多くの拍動数を必要とします）。 |

---

## 🛠️ ソフトウェア構成とモジュール

コードベースはモジュール式のパイプラインとして構築されています：

1.  📂 **[io_loader.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/io_loader.py):** I/O マネージャー。CSV/TXT (Pandas)、MAT (SciPy)、EDF/BDF (pyedflib) ファイルをディスクまたはメモリバッファから読み込み、統一された `SignalData` オブジェクトを返します。
2.  𝄠 **[preprocessing.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/preprocessing.py):** ゼロ位相バターワースバンドパスフィルタ（$0.5 - 30 \text{ Hz}$）、線形デトレンド、拍動セグメント化（BioSPPy）、最小最大振幅スケーリング、128ポイント線形補間。
3.  🌀 **[adaptive_engine.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/adaptive_engine.py):** 中核エンジン。$H_{morph}$ トレンドを計算し、窓幅とオーバーラップを写像し、$T_{crit}$ 制約を適用して拍動インデックスのスライス区間を生成します。
4.  📊 **[feature_extraction.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/feature_extraction.py):** 21個の形態的特徴量（時間、面積、VPG、SDPTGの $a,b,c,d,e$ 点、および加齢指数 AGI）を抽出します。可変窓と標準の30秒固定窓のベンチマークを比較します。
5.  📈 **[visualization.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/visualization.py):** Plotly エンジン。ズームが同期するサブプロットと、Streamlit のデータエクスポート（CSV、JSON）パネルを生成します。
6.  🚀 **[app.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/app.py):** Streamlit インターフェース。実行ステージを統括し、リアルタイム応答のための**段階的キャッシュ（Cascade Caching, `@st.cache_data`）**を管理します。

---

## 🚀 はじめに

### 動作環境
本プロジェクトには Python 3.10+ と以下のパッケージが必要です：
*   `numpy`
*   `pandas`
*   `scipy`
*   `pyedflib`
*   `biosppy`
*   `peakutils`
*   `plotly`
*   `streamlit`

### 1. リポジトリのクローンと仮想環境の設定
```bash
git clone https://github.com/your-username/adaptive-ppg-windows.git
cd adaptive-ppg-windows

# 仮想環境の作成
python -m venv venv
source venv/bin/activate  # Windowsの場合: venv\Scripts\activate
```

### 2. 依存関係のインストール
```bash
pip install numpy pandas scipy pyedflib biosppy peakutils plotly streamlit
```

### 3. ダッシュボードの起動
```bash
streamlit run app.py
```

### 4. モジュール単体テスト of 実行
各ファイルには擬似データを生成するテストブロックが含まれています。以下のコマンドを実行することで単独でテストできます：
```bash
python io_loader.py
python preprocessing.py
python adaptive_engine.py
python feature_extraction.py
python visualization.py
```

---

## 📊 サンプルデータセット

ダッシュボードをテストしやすくするため、`sample_data/` フォルダに臨床記録のサンプルが含まれています。
- **ファイル**: `sample_data/bidmc_01_Signals.csv`
- **信号列 (Channel)**: `PLETH` (PPG 光電脈波を指します)
- **サンプリング周波数**: $125 \text{ Hz}$

### データセットの出典および引用
このサンプルデータは、**PhysioNet** で公開されている **BIDMC PPG and Respiration Dataset** から抽出されたものです：
> *Pimentel, M. A. F., Charlton, P. H., & Clifton, D. A. (2016). BIDMC PPG and Respiration Dataset (version 1.0.0). PhysioNet. https://doi.org/10.13026/C2GG69.*
> 
> *Pimentel, M. A. F., et al. "Towards a Robust Estimation of Respiratory Rate from Pulse Oximeters." IEEE Transactions on Biomedical Engineering, 64(8), pp. 1914-1923, 2016.*

