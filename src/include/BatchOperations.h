// Algorithm 4.1: Batch Vector-Matrix Multiplication 
// a: cleartext vector (bits), B: packed ciphertexts
TensorRGSW VecMatMult(const std::vector<int>& a, const std::vector<TensorRGSW>& B, const BatchParams& bp);

// Algorithm 5.2: Batch-Permute 
// 置換 pi に基づいてスロットを入れ替える
TensorRGSW Batch_Permute(const TensorRGSW& C, const std::vector<int>& pi);

// Algorithm 5.3: Inv-Auto 
// Automorphism map: xi -> xi^-1
TensorRGSW Inv_Auto(const TensorRGSW& C);

// Algorithm 5.4: Batch-Anti-Rot 
// 重要な構成要素: Permute と Inv-Auto を組み合わせて回転を実現
TensorRGSW Batch_Anti_Rot(const TensorRGSW& C, int delta);
