# N05-N06 学習ログ

## N05：抽出前の予想とbaseline

- 残す挙動・実際の重複・N04選択：
  - 【要記入：抽出前に自分で予想したこと】
  - 実際の重複：N01〜N04の各 `run.jl` に、周期隣接の添字、線形流束・Burgers流束、新旧バッファの検証、保存区間に合わせた刻み（`fit_timestep`）がそれぞれ書かれていた。
  - N04の選択モデルは `selected_model = "linear"`（`results/N05/baseline.toml` に記録）。
- baselineを取得したcommitと対象コード：
  - `N05.jl baseline` は共通化前の commit `bcd571e`（N02完了時点）で取得した（baseline.toml の `git_head`）。
  - 対象コードは N01〜N04 の `run.jl`・`provided_support.jl`・`tests.jl` と `src/` の各ファイル（baseline.toml の `[code_sha256]`）。
  - baseline の SHA-256：`cfcce52c5521b5eb94be022507bd92c4f0b3ac446f18a8dd5b92cb25c3c15162`

## N05：採用したdiffと回帰

- 共通APIへ移した処理と、課題に残した処理：
  - commit `163f4f6` で、`periodic_left_index`／`periodic_right_index`、`linear_flux`／`burgers_flux`、`validate_buffers`、`fit_timestep` を `src/ThermofluidExercise.jl` の共通APIにした。N01〜N04の更新関数（`upwind_step!`、`centered_step!`、`nonlinear_upwind_step!`、`diffusion_step!`、`advection_diffusion_step!` など）も同じファイルへ移し、各 `run.jl` の旧入口は共通APIに処理を任せるだけにした。
  - 課題側に残したもの：初期条件、実験条件、診断、作図、N04のモデル選択。
  - 【要記入：上の分け方を選んだ理由】
- テスト・代表結果の差とその理由：
  - `N05.jl verify`（`run.jl` に含まれる）の `regression.toml` では `[differences]` がすべて `0.0` だった。`baseline_sha256` も baseline と一致した。
  - 最終の `git_head` は `57aaab4`。N06実装後の `src/` を含む `code_sha256` で照合した。
  - 差が 0 なのは、更新式と添字・境界規約を変えずに置き場所だけを移したためである。ただし、回帰がないことは「抽出前と同じ」ことを示すだけで、物理的に正しいことは示さない。正しさは各課題のテスト（解析解・保存量）で別に確認している。
- 復帰点を使った場合の採用差分と説明：
  - 【要記入：復帰支援を使ったかどうか。使っていなければ「使用なし」】

## N06：二次元更新と数値検証

- 軸・周期端・角・方向別寄与・合成CFL：
  - `u[i,j]` は `x[i]`、`y[j]` に対応する。1番目の添字が x、2番目の添字が y。
  - 更新式（`src/N06Advection.jl`）：
    `u_new[i,j] = u_old[i,j] - dt/dx*(F(u_old[i,j]) - F(u_old[im,j])) - dt/dy*(G(u_old[i,j]) - G(u_old[i,jm]))`
    - `F = linear_flux(・, cx)`、`G = linear_flux(・, cy)`
    - `im = periodic_left_index(i, nx)`、`jm = periodic_left_index(j, ny)`
  - 周期端：`i=1` の左隣は `nx`、`j=1` の左隣は `ny`。
  - 角：角の点 `(1,1)` は x 方向に `(nx,1)`、y 方向に `(1,ny)` を参照する。x と y を独立に折り返すので、角でも特別な処理はいらない。
  - 全点を `u_old` だけから計算する。`u_new` と `u_old` は別の配列で、`simulate.jl` では1ステップごとに `a, b = b, a` で入れ替える。
  - 合成CFL：`Cx + Cy = cx*dt/dx + cy*dt/dy <= 1` のとき、新しい値は旧い3点 `(1-Cx-Cy)u[i,j] + Cx u[im,j] + Cy u[i,jm]` の凸結合になる。そのため `stable_timestep` は `dt = safety / (cx/dx + cy/dy)` とした。
    - 配布テスト：`(1, 0.75, 0.5, 0.25)` → `0.8/(2+3) = 0.16`
    - 片方向の速度が 0 の場合は、もう一方向だけの1次元CFLになる。
- 複数ステップの保存量、自作テスト：
  - 自作1「複数ステップの保存量」
    - 入力：7×5（非正方）の非対称な場 `sin(0.9i) + 0.3j^2 + 0.1ij`。`cx=1.3`、`cy=0.4`、`dx=0.3`、`dy=0.2`、`dt=stable_timestep(...)` で20ステップ進める。
    - 判定1：総和が保存される（`rtol=1e-13`）。値は O(1)、点数は35、反復は20回なので、丸め誤差は 1e-15 程度で、1e-13 なら十分な余裕がある。
    - 判定2：最小・最大が初期値の範囲を出ない（±1e-13）。
    - 判定3：旧配列が変わっていない。
    - 判定4：`Cx=1, Cy=0` の1ステップが `circshift(u0,(1,0))` と、`Cx=0, Cy=1` の1ステップが `circshift(u0,(0,1))` と一致する（`atol=1e-13`）。
    - 確かめたこと：y 方向の左隣をわざと右隣に変えると、総和と最小・最大のテストは通ってしまい、シフトのテストだけが失敗した。保存量だけでは隣接の向きの誤りは見つからないので、シフトのテストを入れた。
  - 自作2「3格子の解析解誤差と収束」：公式3格子を `simulate_case` で t=1 まで計算した。`exact_field` とのL2誤差について、(a) 正であること、(b) 格子を細かくすると減ること、(c) 観測次数が `0.8 <= p <= 1.2` に入ることを確認した。
  - 公式出力の質量の相対変化（t=0→1）は、3格子とも `0.0` だった。
- 3格子の誤差・観測次数と考察：

  | 格子 | dt | ステップ数/区間 | L2誤差(t=1) | L∞誤差(t=1) | 観測次数 |
  |---|---|---|---|---|---|
  | 40×30 | 0.022727 | 11 | 0.045867 | 0.095696 | — |
  | 80×60 | 0.011364 | 22 | 0.024065 | 0.050010 | 0.931 |
  | 160×120 | 0.005682 | 44 | 0.012333 | 0.025576 | 0.964 |

  - 一次風上法の打切り誤差は O(Δx, Δy, Δt) なので、格子幅を半分にすると誤差もおよそ半分になる（次数 ≈ 1）。格子を細かくするほど次数は 1 に近づいている。
  - 主な誤差は数値拡散である。`fields.png` の差の図では、ピーク（x≈1.5, y≈0.5）で数値解が低く、谷で高い。振幅がならされている。
  - 【要記入：自分の考察】

## N06：保存場の再解析と再作図

- 分散を追加した結果、HDF5の変更前後のSHA-256：
  - `analyze.jl` を単独で再実行する前後とも、`fields.h5` の SHA-256 は `82598bb1581a53ef91317aa85fb5cde6a5c67424931f9e03c33a60785ff24d32` で変わらなかった。`diagnostics_complete = true`。
  - 空間分散（1/N の母分散）：

    | 格子 | t=0 | t=1 |
    |---|---|---|
    | 40×30 | 0.07 | 0.04830 |
    | 80×60 | 0.07 | 0.05811 |
    | 160×120 | 0.07 | 0.06377 |

  - 解析解は形を保ったまま平行移動するだけなので、分散は 0.07 のまま一定である。0.07 は (0.2²+0.3²+0.1²)/2 で求まる。
  - 数値解の分散が減るのは数値拡散のためである。格子が粗いほど大きく減る。
- 表示設定の変更、入力HDF5・summaryのSHA-256：
  - `plot.jl` で `COLORMAP` を `:viridis` から `:plasma` に、`DISPLAY_CASE` を `"n080x060"` から `"n040x030"` に変えた。`COLOR_LIMITS = (0.4, 1.6)` は変えていない。初期場の値域がちょうど 0.4〜1.6 で、初期場・数値解・解析解の3つの図に共通の範囲である。
  - 粗い格子を表示に選んだのは、数値拡散による振幅の減衰が見やすいからである。
  - 再作図の前後とも、入力HDF5 は `82598bb1…ff24d32`、summary.toml は `1f6bcc49e7cb9dccad87090c4137b7dcc84bf9aa9edb4534dbb05cf9eb2e1f4d` で変わらなかった。`plots.toml` の `source_fields_sha256`／`source_summary_sha256` とも一致した。
  - 再計算・再解析・再作図の使い分け：速度や格子数を変えると `simulate.jl` から再計算が必要である。分散のような解析項目の追加は、保存場から `analyze.jl` で再解析すればよい。表示だけの変更なら `plot.jl` で再作図すればよい。
- 保存されなかった時刻からは求められない量：
  - 保存したのは t = 0, 0.25, 0.5, 0.75, 1 の5時刻だけである。その間の時刻の場、最小・最大、分散は保存場から直接は求められない。たとえば分散が保存時刻の間で単調に減るかどうかや、極値をとる時刻は分からない。必要なら保存時刻を増やして再計算する。
  - 時間更新のバッファ `a`, `b` は使い回されている。保存時には `U[:,:,k] = a` で値をコピーしており、参照だけを保存すると後のステップで上書きされてしまう。

## 理解度チェック・提出前確認

- N05・N06それぞれの理解度チェックと振り返り：
  - 【要記入：`understanding-check-N05.txt`／`understanding-check-N06.txt` を LETUS に提出したかどうかと、振り返り（回答は自分で書く）】
- 実行コマンド、テスト結果、成果物、未解決事項：
  - `julia --project=. exercises/N05-N06_common_package_2d_advection/run.jl`：N05の回帰確認とN06の全段階が成功した。
  - `julia --project=. exercises/N05-N06_common_package_2d_advection/analyze.jl` と `plot.jl` を単独で再実行した（上記のハッシュを確認）。
  - `julia --project=. exercises/N05-N06_common_package_2d_advection/tests.jl`：全30件が通過した。内訳は、N05共通API 5、N06配布 9、自作 6+6、保存結果と出自 4。
  - `julia --project=. -e 'using Pkg; Pkg.test()'`：passed（F03-F04 で警告あり。下記参照）。
  - `julia --project=. scripts/course.jl check-results`：サイズの上限内だった。
  - 成果物：`results/N05/baseline.toml`、`results/N05/regression.toml`、`results/N06/` の6ファイル（`fields.h5`、`summary.toml`、`fields.png`、`diagnostics.png`、`convergence.png`、`plots.toml`）。
  - 未解決事項：
    - `course_progress.toml` の `current` が `N02` のままで、`Pkg.test()` は N03 以降を検査していない。N02 の PR をマージした後、進捗を進める必要がある。
    - `Pkg.test()` で F03-F04 の配布テストがエラーになる（完了済みの課題なので警告扱い）。
- AI利用（Claude Code）：
  - 依頼した内容：
    - リポジトリの状況確認
    - N06の `stable_timestep`／`advection_step!` の実装
    - 自作テスト2件の作成
    - `spatial_variance` の実装
    - 公式出力の生成と、表示設定変更後の再作図
    - このログの下書き
  - 採用した提案：
    - 更新式を共通API（`periodic_left_index`、`linear_flux`）で組み立てる形
    - 分散の定義を 1/N とすること。参照値の t=0 の値 0.07 が、振幅から求めた (0.2²+0.3²+0.1²)/2 と一致することで確かめた。
  - 確認した方法：
    - 配布テストの期待値を手計算と照合した（`new[1,1] = 1 - 0.2(1-9) - 0.3(1-3) = 3.2`）。
    - 実装をわざと壊し、自作テストが失敗することを確かめた。
    - `check_complete` で参照値と照合した。
  - 表示設定は、AIが示した3案から自分で選んだ（粗い格子＋`:plasma`）。
  - 【要記入：自分で理解・修正した点、採用しなかった提案があればその理由】
