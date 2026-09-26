# Anisotropy Time-Series Analyzer

プレートリーダーで測定した蛍光異方性(mA)などの時系列データ(Excel)を読み込み、
プロット、デトレンド、ピーク/トラフ検出、アクトグラム、ピーク/トラフの確認・手動修正、
周期・位相の回帰分析、サインカーブフィッティングによる周期・振幅・位相の解析を行う
デスクトップアプリです。

## 入力ファイルの形式

| 列 | 内容 |
|---|---|
| 1列目 | 時間 (h) |
| 2列目 | 温度 (°C) |
| 3列目以降 | 各ウェルの測定値(1行目にウェル名: D4, D5, …) |

行数(測定時間)と列数(ウェル数)はファイルごとに異なってかまいません。

## Windows 用 EXE の入手

### GitHub Actions で自動ビルド

1. **Actions** タブ → **Build Windows EXE** → 最新の実行(緑のチェック)を開く
2. ページ下部の **Artifacts** から **AnisotropyAnalyzer-windows** をダウンロード
3. ZIP を展開すると `AnisotropyAnalyzer.exe` が入っています

ビルドは `main` ブランチへの push のたびに自動で行われます。
Actions タブの **Run workflow** ボタンから手動でも実行できます。

### リリースとして公開する

`v` で始まるタグを push すると、EXE を添付した Release が自動で作成されます。

```
git tag v1.0.0
git push origin v1.0.0
```

(GitHub の Web 画面で **Releases → Draft a new release** からタグ `v1.0.0` を作っても同じです)

### 自分の PC でビルドする

Python 3.10 以上をインストールした Windows で `build_exe_local.bat` をダブルクリックすると、
`dist\AnisotropyAnalyzer.exe` が作られます。

## EXE の使い方と注意

- ダブルクリックで起動します(インストール不要)。
- 起動には 10〜20 秒ほどかかることがあります(1ファイルにまとめた EXE は、起動時に一時フォルダへ展開するため)。初回はさらに時間がかかります。
- 署名のない EXE のため、初回起動時に **「Windows によって PC が保護されました」** と表示されることがあります。**詳細情報 → 実行** で起動できます。
- ウイルス対策ソフトが誤検知することがあります(PyInstaller で作った EXE ではよくある現象です)。

### コマンドラインオプション

| オプション | 内容 |
|---|---|
| `--selftest [フォルダ]` | 合成データで解析と Excel/PDF 出力を実行し、動作を確認します。結果は `フォルダ\selftest.log` に書かれ、終了コード 0 なら成功です。 |
| `--version` | バージョンを表示します。 |

## Python から直接実行する場合

```
pip install -r requirements.txt
python anisotropy_analyzer_app.py
```

## ファイル構成

| ファイル | 内容 |
|---|---|
| `anisotropy_analyzer_app.py` | アプリ本体(1ファイル) |
| `requirements.txt` | 必要なライブラリ |
| `requirements-build.txt` | EXE の作成に必要なライブラリ (PyInstaller) |
| `.github/workflows/build-windows-exe.yml` | GitHub Actions で EXE を作るワークフロー |
| `build_exe_local.bat` | 自分の PC で EXE を作るバッチファイル |
