# F01 学習ログ

## 予想と結果

- 変更した箇所：exercises/F01_first_pull_request/run.jl の student_greeting 内の仮実装（TODO行）を `"Hello, $(normalized_name)!"` に置き換えた。前後空白の除去と空文字拒否は配布済みコードをそのまま利用した。
- 期待した表示と実際の表示：`julia --project=. exercises/F01_first_pull_request/run.jl "Genki Nagaya"` で `Hello, Genki Nagaya!` を期待し、同じ出力を得た。引数なしでは既定値により `Hello, student!` と表示される。

## 自分のテスト

- 自分で選んだ入力と期待値：`"  Euler  "` → `"Hello, Euler!"`
- 保証すること・保証しないこと：保証するのは、前後の半角スペースが除去されること、出力が `Hello, ` と `!` で挟まれた形式であること。保証しないのは、タブ・改行・全角スペースの扱い、空文字列 `""` の拒否、非ASCII文字を含む名前での挙動、main 関数の出力、および ArgumentError のメッセージ文面（`@test_throws` は例外の型のみを検査する）。
- 実行コマンドと結果：`julia --project=. exercises/F01_first_pull_request/tests.jl` で必須2件・自作1件が成功。`julia --project=. -e 'using Pkg; Pkg.test()'` で全テスト成功（test/tests.jl 199件、F01 tests.jl 3件）。

## つまずきと判断

- つまずき、その対処と理由：文字列補間を `"Hello, $normalized_name!"` と書いたところ `UndefVarError: normalized_name! not defined` となった。Julia では `!` が識別子に使える文字であり、`$normalized_name!` が変数 `normalized_name!` として解釈されたことが原因。`$(normalized_name)` と括弧で変数名の範囲を明示して解決した。
- 配布済み必須テストが保証することと、自作テストの入力・期待値を選んだ理由：必須テストは「前後の空白が除去されること」と「空白のみの名前が ArgumentError で拒否されること」の2点のみを検証している。自作テストは必須テストと異なる名前を用いたが、検証している性質は前者と同一であり、新たな保証を追加するものではない。

## AI利用・出典

- 依頼内容：課題ページの読解と手順整理、student_greeting の実装方針、自作テストケースの選定、本学習ログの下書き作成。
- 採用・修正・却下と理由：実装案は採用したが、提案された補間表記 `"Hello, $normalized_name!"` はエラーとなったため、エラーメッセージから原因を特定し `$(normalized_name)` に修正して採用した。
- 外部素材の出典、変更点、利用条件：なし
- リポジトリ作り直しに伴う再提出（2026-11-06）：リポジトリを削除・再作成したため、以前の自分の F01 実装をバックアップから戻して再提出した。この作業には Claude Code を使用し、ファイルの復元、テストと実行の再確認、本ログの修正を依頼した。その際、自作テストが `"Hello, Euler!"` という文字列を置いただけで `@test` による判定をしていないという指摘を受け、`@test F01FirstPullRequest.student_greeting("  Euler  ") == "Hello, Euler!"` に修正して採用した（修正前はテストとして何も検証していなかった）。

## 理解度チェック・LETUS提出

- 対応する授業ID：F01
- LETUS提出日：
- 提出済み確認：
- 理解できた点：git commandの使い方
