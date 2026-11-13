# F03-F04の課題テスト。提供済みの必須確認 → 受講生が追加する自作テストの順に読みます。
# 必須テストの期待値は保持し、自作テストのTODOを編集します。

using Test

if !isdefined(Main, :F04NumericalDifferentiation)
    include(joinpath(@__DIR__, "run.jl"))
end
using .F04NumericalDifferentiation

# 提供済みの期待値を使い、実装と独立した基準で確認する。
@testset "F03-F04 必須テスト（配布済み）" begin
    @testset "差分値（forward_difference・backward_difference・centered_difference）" begin

        """
            quadratic(x)

        必須テストで使う2次関数を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        quadratic(x) = x ^ 2
        @test forward_difference(quadratic, 2.0, 0.5) ≈ 4.5
        @test backward_difference(quadratic, 2.0, 0.5) ≈ 3.5
        @test centered_difference(quadratic, 2.0, 0.5) ≈ 4.0
    end

    @testset "収束次数（convergence_study）" begin

        """
            cubic(x)

        必須テストで使う3次関数を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        cubic(x) = x ^ 3

        """
            cubic_derivative(x)

        必須テストで使う3次関数の解析的な微分を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        cubic_derivative(x) = 3x ^ 2
        study = convergence_study(cubic, cubic_derivative, 1.0, [0.2, 0.1, 0.05])
        @test all(ratio -> 1.8 <= ratio <= 2.2, study.forward_ratios)
        @test all(ratio -> 1.8 <= ratio <= 2.2, study.backward_ratios)
        @test all(ratio -> 3.9 <= ratio <= 4.1, study.centered_ratios)
    end

    @testset "ベクトル公式の格子収束（verify_vector_identities）" begin
        coarse = verify_vector_identities(9)
        fine = verify_vector_identities(17)
        @test keys(coarse) ==
              (:curl_gradient, :divergence_curl, :product_divergence, :curl_curl)
        for key in keys(coarse)
            @test 0 < fine[key] < coarse[key]
            @test 3.0 <= coarse[key] / fine[key] <= 4.8
        end
    end
end

# 必須テストが扱わない条件を選び、期待値を自分で決める。
@testset "F03-F04 自作テスト" begin
    # TODO(自作): 別の関数、評価点、または入力条件を選び、どの実装ミスを検出するか説明できるテストを一つ書く。F03単独の別テストは作らない。
    # 入力: f(x)=x^3、x=1、h=0.5（必須テストの2次関数では中心差分の誤差が0になり、誤差項を確かめられない）。
    # 期待値（手計算）: 前進 (1.5^3-1)/0.5=4.75、後退 (1-0.5^3)/0.5=1.75、中心 (1.5^3-0.5^3)/(2*0.5)=3.25。
    # 3次式ではTaylor展開が有限で打ち切られ、中心差分の誤差は h^2 f'''/6 = 0.25*6/6 = 0.25 にちょうど等しい。
    # 検出する誤り: 中心差分の分母をhにする(6.5)、前進と後退の取り違え(4.75と1.75が入れ替わる)、
    # 符号の反転、中心差分を前進と後退の平均ではなく片側で計算する誤り。
    # 許容値: 値はO(1)で演算は数回なので、rtol=1e-14で丸め誤差だけを許す。
    cubic3(x) = x ^ 3
    @test forward_difference(cubic3, 1.0, 0.5) ≈ 4.75 rtol = 1e-14
    @test backward_difference(cubic3, 1.0, 0.5) ≈ 1.75 rtol = 1e-14
    @test centered_difference(cubic3, 1.0, 0.5) ≈ 3.25 rtol = 1e-14
    @test centered_difference(cubic3, 1.0, 0.5) - 3.0 ≈ 0.5 ^ 2 * 6 / 6 rtol = 1e-14
    # 中心差分は前進と後退の平均に等しい（3つの実装の整合）。
    @test centered_difference(cubic3, 1.0, 0.5) ≈
          (forward_difference(cubic3, 1.0, 0.5) + backward_difference(cubic3, 1.0, 0.5)) / 2 rtol =
        1e-14
end
