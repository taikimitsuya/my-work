class TensorRGSW {
public:
    // 内部的には MK-TFHE の TRGSW や vector<TLwe> を保持
    std::vector<TRGSW> ciphertexts;
    BatchParams::Mode current_mode; 

    // Constructor
    TensorRGSW(const BatchParams& params);

    // Algorithm: RGSW-Pack 
    // inputs: ベクトル化されたRGSW暗号文
    void RGSW_Pack(const std::vector<TRGSW>& inputs, BatchParams::Mode target_mode);

    // Algorithm: UnPack 
    std::vector<TRGSW> UnPack();
    
    // Homomorphic Matrix-Vector Multiplication building block
    // Algorithm: Batch-Mult 
    // モード間の遷移ルール (R12 * (R12->R13) => R13 など) を実装
    friend TensorRGSW Batch_Mult(const TensorRGSW& A, const TensorRGSW& B);
};
