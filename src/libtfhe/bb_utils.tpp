// bb_utils.tpp

namespace bbii {

template <typename T>
std::vector<T> rearrange(const std::vector<T>& input, int32_t d) {
    // 論文 [cite: 627] "Rearr: on input m -> 2d ... outputs (coeffs(a'_1)...)"
    // これは多項式の係数を、再帰的な構造に合わせて並べ替える処理。
    // 一般的な解釈として、係数ベクトルをストライドアクセスで再配置する。
    // 例: a(X) = a0 + a1*X + ... -> a(Y, Z) への変換
    
    // 実装: ストライド置換 (Stride Permutation)
    // 入力サイズ N を d 個のブロックに分けるイメージ
    
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be divisible by d");
    
    std::vector<T> output(n);
    size_t rows = d;
    size_t cols = n / d;
    
    // ストライド順序で配置
    // input[i * cols + j] -> output[j * rows + i] 
    // これは行列の転置(Transpose)のような操作
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            output[j * rows + i] = input[i * cols + j];
        }
    }
    
    return output;
}

template <typename T>
std::vector<T> reverse_rearrange(const std::vector<T>& input, int32_t d) {
    // 論文 [cite: 628] "Rev-Rearr: on input 2d -> m"
    // Rearr の逆操作
    
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be divisible by d");

    std::vector<T> output(n);
    size_t rows = d; // 元の rows
    size_t cols = n / d; // 元の cols
    
    // 逆転置
    // input[j * rows + i] -> output[i * cols + j]
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            output[i * cols + j] = input[j * rows + i];
        }
    }
    
    return output;
}

} // namespace bbii