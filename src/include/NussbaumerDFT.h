// Algorithm 5.5: RGSW.EncVec-MatMult 
// 暗号化されたベクトルと、1の冪乗根を持つ特殊行列 M との乗算
// 実際には DFT行列や逆DFT行列との積に使用される
std::vector<TensorRGSW> EncVec_MatMult(
    const Matrix& M, // Entries are powers of xi
    const std::vector<TensorRGSW>& C_in,
    const BatchParams& bp
);

// Algorithm 6.3: Hom-DFT^-1 (Recursive) 
// 再帰的にNussbaumer変換を適用する
std::vector<TensorRGSW> Hom_DFT_Inverse(
    int rho_prime, 
    int n_prime, 
    const std::vector<TensorRGSW>& C_in,
    const BatchParams& bp
);
