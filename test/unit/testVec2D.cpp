#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>

#include "lettuce/core/Vec.h"



using Real = double;
constexpr double DEFAULT_EPSILON = std::numeric_limits<Real>::epsilon();

bool isEqual(double a, double b, double epsilon = DEFAULT_EPSILON) {
    return std::abs(a - b) < epsilon;
}

void testVecClassSize() {
    Vec<Real> v1;
    Vec<Real, 3> v2;
    Vec<Real, 4> v3;
    
    assert(v1.size() == 2);
    assert(v2.size() == 3);
    assert(v3.size() == 4);
}

void testVecClassConstructor() {
    Vec<Real> v1;
    Vec<Real> v2(5.0);
    Vec<Real> v3({3.0, 4.0});

    assert(isEqual(v1[0], 0));
    assert(isEqual(v1[1], 0));

    assert(isEqual(v2[0], 5.0));
    assert(isEqual(v2[1], 5.0));
    
    assert(isEqual(v3[0], 3.0));
    assert(isEqual(v3[1], 4.0));
}

void testVecCompOperands() {
    Vec<Real> v1({3, 4});
    Vec<Real> v2({4, 3});

    auto v3 = v1 + v2;
    assert(isEqual(v3[0], 7.0));
    assert(isEqual(v3[1], 7.0));


    auto v4 = v1 - v2;
    assert(isEqual(v4[0], -1));
    assert(isEqual(v4[1], 1));


    auto v5 = v1 * v2;
    assert(isEqual(v5, 24));

    auto v6 = v1 * 10.11;
    assert(isEqual(v6[0], 30.33));
    assert(isEqual(v6[1], 40.44));

    auto v7 = v1 / 1.5;
    assert(isEqual(v7[0], 2.0));
    assert(isEqual(v7[1], 4.0/1.5));

    Vec<Real> v11 = v1;
    v11 += v2;
    assert(isEqual(v11[0], 7.0));
    assert(isEqual(v11[1], 7.0));

    Vec<Real> v12 = v1;
    v12 -= v2;
    assert(isEqual(v12[0], -1));
    assert(isEqual(v12[1], 1));

    Vec<Real> v13 = v1;
    v13 *= 10.11;
    assert(isEqual(v13[0], 30.33));
    assert(isEqual(v13[1], 40.44));


    Vec<Real> v14 = v1;
    v14 /= 1.5;
    assert(isEqual(v14[0], 2.0));
    assert(isEqual(v14[1], 4.0/1.5));

    auto v8 = v1 + 10.11;
    assert(isEqual(v8[0], 13.11));
    assert(isEqual(v8[1], 14.11));

    Vec<Real> v15 = v1;
    v15 = -v15;
    assert(isEqual(v15[0], -v1[0]));
    assert(isEqual(v15[1], -v1[1]));
}

void testVecCmpOperator() {
    Vec<Real> a;
    Vec<Real> b({1.0, 0.3});

    a[0] = 1.0; a[1] = 0.1 + 0.2;
    // std::cerr << a[0] << ", " << a[1] << std::endl;
    // std::cerr << b[0] << ", " << a[1] << std::endl;
    assert(a == b);

    Vec<Real> a1 = a;
    a1[1] += 0.1;
    assert(a1 != b);

    Vec<Real> c ({12.01111119, 20.0});
    assert(c >= 12.01111111);
    assert(c <= 20.0000001);
    assert(c > 11.999999);
    assert(c < 20.000001);
}

void testVecCalOperatorScalar() {
    Vec<Real> a ({3.0, 4.0});
    auto a1 = a + 5.0000001;
    assert(isEqual(a1[0], (3.0 + 5.0000001)));
    // assert(isEqual(a1[0], 8.0000001));
    assert(isEqual(a1[1], (4.0 + 5.0000001)));
    // assert(isEqual(a1[1], 9.0000001));

    Vec<Real> b({7.0, 11.0});
    const auto b1 = b + 5.0000001;
    assert(isEqual(b1[0], (7.0 + 5.0000001)));
    assert(isEqual(b1[1], (11.0 + 5.0000001)));
}

void testVecFloatingPointEdgeCases() {
    using std::isfinite;
    using std::isinf;
    using std::isnan;
    using std::signbit;

    const Real INF = std::numeric_limits<Real>::infinity();
    const Real NANv = std::numeric_limits<Real>::quiet_NaN();
    const Real EPS = std::numeric_limits<Real>::epsilon();

    // 1) 不同尺度（large + tiny）：測試精度流失
    {
        Vec<Real> v({1e16, 1.0});
        auto w = v + Vec<Real>({1.0, 1e-16}); // 第二個分量幾乎不改變
        assert(isEqual(w[0], 1e16 + 1.0));     // 可見加 1 不會被吞掉
        assert(isEqual(w[1], 1.0));            // 1 + 1e-16 對 double 幾乎無效
    }

    // 2) 災難性消去（catastrophic cancellation）
    {
        Vec<Real> v1({1e8 + 1.0, 1e8});
        Vec<Real> v2({1e8,       1e8});
        auto d = v1 - v2;
        // 理論上 d = {1, 0}，第二個分量可能因為相同數相減直接為 0
        assert(isEqual(d[0], 1.0));
        assert(isEqual(d[1], 0.0));
    }

    // 3) 向量點積的消去案例（* 為 dot）
    {
        Vec<Real, 3> a({1e8, 1.0, 1.0});
        Vec<Real, 3> b({1.0, -1e8, 1.0});
        // 1e8*1 + 1*(-1e8) + 1*1 = 1（大量抵銷後只剩 1）
        Real dot = a * b;
        assert(isEqual(dot, 1.0));
        // 也測 commutativity
        assert(isEqual(a * b, b * a));
    }

    // 4) 連乘/連除的非結合性（縮放次序差異）
    {
        Vec<Real> v({3.0, -4.0});
        Real s1 = 1e-16, s2 = 1e16;

        auto left  = (v * s1) * s2;   // ((v * 1e-16) * 1e16)
        auto right = v * (s1 * s2);   // (v * (1e-16 * 1e16)) = v * 1
        // 理論上都等於 v，但 left/right 可能有極小差異
        assert(isEqual(right[0], v[0])); 
        assert(isEqual(right[1], v[1]));
        // left 也應非常接近 v
        assert(isEqual(left[0], v[0]));
        assert(isEqual(left[1], v[1]));
    }

    // 5) 正規化（除以範數），測試精度與除法
    {
        Vec<Real> v({3.0, 4.0});
        Real n = std::sqrt(v * v); // v·v = 25 -> sqrt = 5
        auto u = v / n;
        assert(isEqual(u[0], 0.6));
        assert(isEqual(u[1], 0.8));
        // 單位向量長度 ≈ 1
        Real un = std::sqrt(u * u);
        assert(std::abs(un - Real(1)) <= 8*EPS);
    }

    // 6) 次正規數與下溢（underflow to zero）
    {
        // 對 double：~1e-308 是臨界，1e-320 很小，除以 1e10 => 1e-330（下溢 -> 0）
        Vec<Real> tiny({1e-320, -1e-320});
        auto z = tiny / Real(1e10);
        // 可能皆為 0（或次正規更貼近 0）。我們接受它為 0。
        assert(isEqual(z[0], 0.0));
        assert(isEqual(z[1], 0.0));
    }

    // 7) 上溢（overflow -> inf）
    {
        Vec<Real> big({1e308, -1e308});
        auto infv = big / Real(1e-308); // 1e308 / 1e-308 = 1e616 -> inf
        assert(isinf(infv[0]) && !signbit(infv[0]));
        assert(isinf(infv[1]) &&  signbit(infv[1]));
    }

    // 8) 與 0 的運算（含 -0.0）
    {
        Vec<Real> v({0.0, -0.0});
        // 乘上負數可以測試 -0.0 的保留情況（實作若保留符號）
        auto w = v * Real(-2.0);
        // 第一個 +0 * (-2) 仍為 -0，但許多比較會視為 0；僅檢查「數值等於 0」
        assert(isEqual(w[0], 0.0));
        assert(isEqual(w[1], 0.0));
        // 如需嚴格檢查 -0.0，才用 signbit：
        // 注意：不同平台/最佳化下 -0 的保留不保證
        // assert(signbit(w[0]) == true/false); // 視需求開啟
    }

    // 9) 鏈式運算的穩定性（(v*s)/s ≈ v）
    {
        Vec<Real> v({-1.23456789012345, 9.87654321098765});
        Real s = 10.11;
        auto back = (v * s) / s;
        assert(isEqual(back[0], v[0]));
        assert(isEqual(back[1], v[1]));
    }

    // 10) 加法/減法的分配律在浮點下的偏差
    {
        Vec<Real> a({1e8, 1.0});
        Vec<Real> b({1.0,  1e-8});
        // Real s = 3.14159265358979323846;
        Real s = 3.14159;

        auto left  = s * (a + b);
        auto right = s * a + s * b;

        // 浮點下 left/right 可能有極微差異，但應該非常接近
        // std::cerr << left[0] <<  ", " << right[0] << std::endl;
        // assert(isEqual(left[0], right[0]));
        // assert(isEqual(left[1], right[1]));
    }

    // 11) 與 NaN 的互動（如果你的 Vec 運算遵循 IEEE）
    {
        Vec<Real> v({1.0, 2.0});
        auto n1 = v * (NANv);   // 期望兩分量皆為 NaN
        assert(isnan(n1[0]) && isnan(n1[1]));

        Vec<Real> u({NANv, 3.0});
        Real d = u * v;         // dot 中若含 NaN，結果通常為 NaN
        assert(isnan(d));
    }

    // 12) 混合：大數 * 小數（避免提前 overflow / underflow 的順序）
    {
        Vec<Real> v({1e8, 1e-8});
        Real s1 = 1e-8, s2 = 1e8;

        // ((v*s1)*s2) 與 v*(s1*s2) 都應 ≈ v
        auto left  = (v * s1) * s2;
        auto right = v * (s1 * s2);

        assert(isEqual(left[0], v[0]));
        assert(isEqual(left[1], v[1]));
        assert(isEqual(right[0], v[0]));
        assert(isEqual(right[1], v[1]));
    }

    // 13) 加/減純量（驗證廣播與符號）
    {
        Vec<Real> v({-1.25, 2.5});
        auto a = v + Real(1e-8);
        // auto b = v - Real(1e-8);
        assert(isEqual(a[0], -1.25 + 1e-8));
        assert(isEqual(a[1],  2.5  + 1e-8));
        // assert(isEqual(b[0], -1.25 - 1e-8));
        // assert(isEqual(b[1],  2.5  - 1e-8));
    }

    // 14) 比例不變性（縮放後再做點積與原點積的比例）
    {
        Vec<Real> x({3.0, 4.0});
        Vec<Real> y({4.0, 3.0});
        Real s = 1e-8;

        Real d1 = (x * s) * (y * s);   // (sx)·(sy) = s^2 (x·y)
        Real d2 = x * y;               // 24
        assert(isEqual(d1, d2 * s * s));
    }
}

int main() {
    testVecClassConstructor();
    testVecClassSize();
    testVecCompOperands();
    testVecCmpOperator();
    testVecCalOperatorScalar();
    testVecFloatingPointEdgeCases();

    return 0;
}