// SPDX-License-Identifier: MIT
pragma solidity >=0.4.22 <0.9.0;
    // This smart contract estimtes the gas cost of the verification functions (PVADVerification and fairnessVerification). We omit simple functions such as staking, submitting an opening of a polynomical commtiment etc in this version.
contract VerifyPVAD {
    uint256 constant PRIME_Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 constant BABYJUB_P = 21888242871839275222246405745257275088548364400416034343698204186575808495617;

    struct G1Point {
        uint256 X;
        uint256 Y;
    }

    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint256[2] X;
        uint256[2] Y;
    }

    struct G1G2 {
        G1Point g1;
        G2Point g2;
    }

    struct zkPoCProof{
        G1Point g_one_over_xPrime; // g_1/x'
        G2Point h_alpha_rPrime; // h_αr'
        uint256 zx;
        uint256 zr; 
        uint256 beta;
    }

    // h
    G2Point SRS_h = G2Point(
        [11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781],
        [4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930]
    );
    // h^α
    G2Point SRS_hAlpha = G2Point(
        [11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781],
        [4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930]
    );

    // h^beta
    G2Point SRS_hBeta = G2Point(
        [11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781],
        [4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930]
    );
    // h^αx
    G2Point SRS_hAlphaX = G2Point(
        [11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781],
        [4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930]
    );

    // The G1 generator; g
    function SRS_g() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }

    uint256 yz = mulmod(y, z, BABYJUB_P);
    uint256 z = uint256(2);
    uint256 y = uint256(3);
    uint256 zd_prime = uint256(4);
    uint256 beta = uint256(1);

    // number of constraints

    uint256 N1 = 232;
    uint256 N2 = 944;
    uint256 N3 = 3880;
    uint256 N4 = 17336;
    // length of input
    uint256 d_i_plus_one = 12;
    // length of algorithm 
    uint256 d_f_1 = 22; 
    uint256 d_f_2 = 10; 
    uint256 d_f_3 = 22; 
    uint256 d_f_4 = 46; 
    // length of output
    uint256 d_o_plus_one = 2;
    
    // t_2
    uint256 deadline = 11559732032986387107991004021392285783925812861821192530917403151452391805634;

    // Commitment from the setup function
    G1Point K = G1Point(
        uint256(20435686948508171234472206488737953800505595616105823290561271581793730135986),
        uint256(7613038940582986439878577004424311309737615170791456916446723479068371769225)
    );
    
    G1Point S_y = G1Point(
        uint256(20435686948508171234472206488737953800505595616105823290561271581793730135986),
        uint256(7613038940582986439878577004424311309737615170791456916446723479068371769225)
    );

    // Commitment from DM
    G1Point F = G1Point(
        uint256(20435686948508171234472206488737953800505595616105823290561271581793730135986),
        uint256(7613038940582986439878577004424311309737615170791456916446723479068371769225)
    );

    // Commitment from user
    G1Point I = G1Point(
        uint256(20435686948508171234472206488737953800505595616105823290561271581793730135986),
        uint256(7613038940582986439878577004424311309737615170791456916446723479068371769225)
    );

    // pk1 = g^(1/-z')
    G1Point pk1 = G1Point(
        uint256(20435686948508171234472206488737953800505595616105823290561271581793730135986),
        uint256(7613038940582986439878577004424311309737615170791456916446723479068371769225)
    );

    // pk2 = h^alpha^(-z')
    G2Point pk2 = G2Point(
        [11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781],
        [4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930]
    );

    struct DM_Proof{
        zkPoCProof pi_Z_prime;
        G1Point pi_R_tilde; 
        G1Point pi_R1; 
        G1Point pi_R2; 
        G1Point pi_T; 
        G1Point pi_K; 
        G1Point pi_S_hat;
        G1Point pi_S1;
        G1Point pi_S2;
        G1Point pi_I_zeta;
        G1Point pi_F_zeta;
        G1Point pi_O_d_prime;
        G1Point pi_O_zeta;
        G1Point pi_O_z_prime;
        uint256 r_tilde; 
        uint256 r_1; 
        uint256 r_2; 
        uint256 t; 
        uint256 k;
        uint256 s_hat;
        uint256 s_1; 
        uint256 s_2; 
        uint256 i_zeta;
        uint256 f_zeta;
        uint256 o_d_prime;
        uint256 o_zeta;
        uint256 o_z_prime;
        // Commitments
        G1Point R_tilde;
        G1Point R;
        G1Point T;
        G1Point S_x;
        G1Point O;
        G1Point O_prime;    
    }

    struct specialProof{
        G1Point _commitment; // F
        G1Point _proof; // pi
        G2Point h_negatealphaz; // h_-αz
        G1Point _index;  // g_nagate1/z
        uint256 _value;  // v
    }

    struct CheckDMax{
        G1Point _com1; // I
        G1Point _com2; // I'
        G2Point _h1; // h^beta
        G2Point _h2; // h^ax(di-d) 
    }

    // rKZG single-point evaluation
    function KZGVerify(
        G1Point memory _commitment, // F
        G1Point memory _proof, // W or pi
        uint256 _index,  // z
        uint256 _value  // F(z) or v
    ) public view returns (bool) {
        // Make sure each parameter is less than the prime q
        require(_commitment.X < BABYJUB_P, "Verifier.verifyKZG: _commitment.X is out of range");
        require(_commitment.Y < BABYJUB_P, "Verifier.verifyKZG: _commitment.Y is out of range");
        require(_proof.X < BABYJUB_P, "Verifier.verifyKZG: _proof.X is out of range");
        require(_proof.Y < BABYJUB_P, "Verifier.verifyKZG: _proof.Y is out of range");
        require(_index < BABYJUB_P, "Verifier.verifyKZG: _index is out of range");
        require(_value < BABYJUB_P, "Verifier.verifyKZG: _value is out of range");
       
        G1Point memory negProof = negate(mulScalar(_proof, _index));
        G1Point memory mulProof = plus(mulScalar(SRS_g(), _value), negProof);

        return pairing_3point(_proof, SRS_hAlphaX,
                                mulProof, SRS_hAlpha,
                                negate(_commitment), SRS_h);
    }

    // rKZG speical verification
    function kzgSpecialVerification(
        specialProof memory cm
    ) public view returns (bool) {
        // Make sure each parameter is less than the prime q
        require(cm._commitment.X < BABYJUB_P, "Verifier.verifyKZG: _commitment.X is out of range");
        require(cm._commitment.Y < BABYJUB_P, "Verifier.verifyKZG: _commitment.Y is out of range");
        require(cm._proof.X < BABYJUB_P, "Verifier.verifyKZG: _proof.X is out of range");
        require(cm._proof.Y < BABYJUB_P, "Verifier.verifyKZG: _proof.Y is out of range");
        require(cm.h_negatealphaz.X[0] < BABYJUB_P, "Verifier.verifyKZG: h_negatealphaz.X[0] is out of range");
        require(cm.h_negatealphaz.Y[0] < BABYJUB_P, "Verifier.verifyKZG: h_negatealphaz.Y[1] is out of range");
        require(cm.h_negatealphaz.X[1] < BABYJUB_P, "Verifier.verifyKZG: h_negatealphaz.X[0] is out of range");
        require(cm.h_negatealphaz.Y[1] < BABYJUB_P, "Verifier.verifyKZG: h_negatealphaz.Y[1] is out of range");
        require(cm._index.X < BABYJUB_P, "Verifier.verifyKZG: g_nagate1/z.X is out of range");
        require(cm._index.Y < BABYJUB_P, "Verifier.verifyKZG: g_nagate1/z.Y is out of range");
        require(cm._value < BABYJUB_P, "Verifier.verifyKZG: _value is out of range");
       
        G1Point memory mulProof = plus(mulScalar(cm._index, cm._value), cm._proof);

        if (pairing_3point(cm._proof, SRS_hAlphaX,
                                mulProof, cm.h_negatealphaz,
                                negate(cm._commitment), SRS_h)){
            return true;
            }
        else return false;
    }

    function zkPoCVerification(
        zkPoCProof memory proof
    ) public view returns (bool) {
        G1Point memory g_zx = mulScalar(SRS_g(), proof.zx); //
        G1Point memory rhs1 = plus(proof.g_one_over_xPrime, mulScalar(pk1, proof.beta));
        G2Point memory h_alphazr = mulScalar2(SRS_hAlpha, proof.zr); //
        G2Point memory rhs2 = add2(proof.h_alpha_rPrime, mulScalar2(SRS_hAlphaX, beta));
        if (g_zx.X == rhs1.X && g_zx.Y == rhs1.Y && h_alphazr.X[1] == rhs2.X[1]&& h_alphazr.Y[0] == rhs2.X[0]&& h_alphazr.Y[0] == rhs2.X[0]&& h_alphazr.Y[1] == rhs2.Y[1]&& proof.beta == uint256(hashFourPoints(SRS_g(),SRS_hAlpha,proof.g_one_over_xPrime,proof.h_alpha_rPrime))){
            return true;
        }
        else return false;
    }
    

    //max degree verification
    function verifyDmax(
        G1Point memory F1,
        G1Point memory F2,
        uint256 d_minus_m
    ) public view returns (bool) {
        if (pairing_2point(F1,SRS_hBeta,F2,mulScalar2(SRS_hAlphaX,d_minus_m))){
            return true;
        }
        else return false;

    }

     // PVAD verification function
    function PVADVerification(
        DM_Proof memory cm
    ) public view returns (bool) {
        bool result = block.timestamp < deadline
        && verifyDmax(cm.O,cm.O_prime,1)
        && zkPoCVerification(cm.pi_Z_prime)
        && KZGVerify(cm.R_tilde,
                      cm.pi_R_tilde,
                      z,
                      cm.r_tilde)
        && KZGVerify(cm.R,
                      cm.pi_R1,
                      z,
                      cm.r_1)
        && KZGVerify(cm.R,
                      cm.pi_R2,
                      yz,
                      cm.r_2)
        && KZGVerify(cm.T,
                      cm.pi_T,
                      z,
                      cm.t)
        && KZGVerify(K,
                      cm.pi_K,
                      y,
                      cm.k)
        && KZGVerify(cm.S_x,
                      cm.pi_S_hat,
                      z,
                      cm.s_hat)
        && KZGVerify(cm.S_x,
                      cm.pi_S1,
                      1,
                      cm.s_1)
        && KZGVerify(S_y,
                      cm.pi_S2,
                      y,
                      cm.s_2)
        && KZGVerify(cm.O,
                      cm.pi_O_d_prime,
                      zd_prime, 
                      cm.o_d_prime)
        && KZGVerify(I,
                      cm.pi_I_zeta,
                      z, 
                      cm.i_zeta)
        && KZGVerify(F,
                      cm.pi_F_zeta,
                      z, 
                      cm.f_zeta) 
        && KZGVerify(cm.O,
                      cm.pi_O_zeta,
                      z, 
                      cm.o_zeta) 
        && y == uint256(hashTwoG1Points(cm.R,cm.R_tilde))
        && cm.t == addmod(mulmod(cm.r_1, 
                                  addmod(cm.r_2,
                                         cm.s_hat, BABYJUB_P), BABYJUB_P),
                            (BABYJUB_P - cm.k), BABYJUB_P)
        && cm.s_1 == cm.s_2
        && z == uint256(hashFourG1Points(cm.R,cm.R_tilde,cm.T,cm.S_x));
        uint256 N_di = addmod(N1, d_i_plus_one, BABYJUB_P);
        uint256 z_n_i = expMod(z, N_di, BABYJUB_P);
        uint256 N_df = addmod(N1, addmod(d_i_plus_one,d_f_1,BABYJUB_P), BABYJUB_P);
        uint256 z_n_f = expMod(z, N_df, BABYJUB_P);
        uint256 N_do = addmod(N1, addmod(d_i_plus_one,addmod(d_f_1,d_o_plus_one,BABYJUB_P),BABYJUB_P), BABYJUB_P);
        uint256 z_n_o = expMod(z, N_do, BABYJUB_P);
        uint256 iPart = mulmod(cm.i_zeta, z_n_i, BABYJUB_P);
        uint256 fPart = mulmod(cm.f_zeta, z_n_f, BABYJUB_P);
        uint256 oPart = mulmod(cm.o_zeta, z_n_o, BABYJUB_P);
        result = result && cm.r_1 == addmod(cm.r_tilde, mulmod(iPart, addmod(fPart,oPart,BABYJUB_P), BABYJUB_P), BABYJUB_P)
        && kzgSpecialVerification(specialProof(cm.O,cm.pi_O_z_prime,pk2,pk1,cm.o_z_prime));
        return result;
    }

    struct fairness_Proof{
        G1Point pi_R_tilde; 
        G1Point pi_R1; 
        G1Point pi_R2; 
        G1Point pi_T; 
        G1Point pi_K; 
        G1Point pi_S_hat;
        G1Point pi_S1;
        G1Point pi_S2;
        G1Point pi_I;
        G1Point pi_O;
        uint256 r_tilde; 
        uint256 r_1; 
        uint256 r_2; 
        uint256 t; 
        uint256 k;
        uint256 s_hat;
        uint256 s_1; 
        uint256 s_2; 
        uint256 i;
        uint256 o;
        // Commitments
        G1Point R_tilde;
        G1Point R;
        G1Point T;
        G1Point S_x;
        G1Point I;
        G1Point O;    
    }

    uint N = 5; // For gas estimation only
    function fairnessVerification(
        fairness_Proof memory cm
        
    ) public view returns (bool) {
        bool  result = KZGVerify(cm.R_tilde,
                      cm.pi_R_tilde,
                      z,
                      cm.r_tilde)
        && KZGVerify(cm.R,
                      cm.pi_R1,
                      z,
                      cm.r_1)
        && KZGVerify(cm.R,
                      cm.pi_R2,
                      yz,
                      cm.r_2)
        && KZGVerify(cm.T,
                      cm.pi_T,
                      z,
                      cm.t)
        && KZGVerify(K,
                      cm.pi_K,
                      y,
                      cm.k)
        && KZGVerify(cm.S_x,
                      cm.pi_S_hat,
                      z,
                      cm.s_hat)
        && KZGVerify(cm.S_x,
                      cm.pi_S1,
                      1,
                      cm.s_1)
        && KZGVerify(S_y,
                      cm.pi_S2,
                      y,
                      cm.s_2)
        
        && cm.t == addmod(mulmod(cm.r_1, 
                                  addmod(cm.r_2,
                                         cm.s_hat, BABYJUB_P), BABYJUB_P),
                            (BABYJUB_P - cm.k), BABYJUB_P)
        && cm.s_1 == cm.s_2;
        // This is for gas estimation only, replace this part with actual input and output commitments
        for (uint i = 0; i < N; i++) {
            result = result && KZGVerify(cm.O,
                      cm.pi_O,
                      z, 
                      cm.o)
            && KZGVerify(I,
                      cm.pi_I,
                      z, 
                      cm.i);
        }
        uint256 N_di = addmod(N1, d_i_plus_one, BABYJUB_P);
        uint256 z_n_i = expMod(z, N_di, BABYJUB_P);
        uint256 iPart;
        uint256 N_do = addmod(N1, addmod(d_i_plus_one,addmod(d_f_1,d_o_plus_one,BABYJUB_P),BABYJUB_P), BABYJUB_P);
        uint256 z_n_o = expMod(z, N_do, BABYJUB_P);
        uint256 oPart;
        for (uint i = 0; i < N; i++) {
            iPart = mulmod(cm.i, z_n_i, BABYJUB_P);
            oPart = mulmod(cm.o, z_n_o, BABYJUB_P);
        }
        
        result = result && cm.r_1 == addmod(cm.r_tilde, addmod(iPart, oPart, BABYJUB_P), BABYJUB_P);
        
        return result;
    }



    //////////////// ********************* //////////////////////////
    // The followings are all Math and pairing functions:

    // Hashes two G1 points by concatenating their X and Y coords
    function hashTwoG1Points(
        G1Point memory p1,
        G1Point memory p2
    ) public pure returns (bytes32) {
        // packs 32+32 bytes for p1, then 32+32 bytes for p2 → 128 bytes total
        return keccak256(
            abi.encodePacked(
                p1.X, p1.Y,
                p2.X, p2.Y
            )
        );
    }
    
    // Hashes four G1 points by concatenating their X and Y coords
    function hashFourG1Points(
        G1Point memory p1,
        G1Point memory p2,
        G1Point memory p3,
        G1Point memory p4
    ) public pure returns (bytes32) {
        // packs 32+32 bytes for p1, then 32+32 bytes for p2 → 128 bytes total
        return keccak256(
            abi.encodePacked(
                p1.X, p1.Y,
                p2.X, p2.Y,
                p3.X, p3.Y,
                p4.X, p4.Y
            )
        );
    }

    function hashFourPoints(
        G1Point memory a,
        G2Point memory b,
        G1Point memory c,
        G2Point memory d
    ) public pure returns (bytes32) {
        return keccak256(
            abi.encodePacked(
                // first G1
                a.X, a.Y,
                // first G2
                b.X[0], b.X[1],
                b.Y[0], b.Y[1],
                // second G1
                c.X, c.Y,
                // second G2
                d.X[0], d.X[1],
                d.Y[0], d.Y[1]
            )
        );
    }

    function ethMessageHash(string memory rawCommitment) internal pure returns (bytes32) {
        return keccak256(
            abi.encodePacked("\x19Ethereum Signed Message:\n32", keccak256(abi.encodePacked(rawCommitment)))
        );
    }
    /*
     * @return The negation of p, i.e. p.plus(p.negate()) should be zero. 
     */
    function negate(G1Point memory p) internal pure returns (G1Point memory) {

        // The prime q in the base field F_q for G1
        if (p.X == 0 && p.Y == 0) {
            return G1Point(0, 0);
        } else {
            return G1Point(p.X, PRIME_Q - (p.Y % PRIME_Q));
        }
    }

    /*
     * @return The sum of two points of G1
     */
    function plus(
        G1Point memory p1,
        G1Point memory p2
    ) internal view returns (G1Point memory r) {

        uint256[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;

        // solium-disable-next-line security/no-inline-assembly
        assembly {
            success := staticcall(sub(gas(), 2000), 6, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }

        require(success, "pairing-add-failed");
    }


    /*
     * @return The product of a point on G1 and a scalar, i.e.
     *         p == p.scalar_mul(1) and p.plus(p) == p.scalar_mul(2) for all
     *         points p.
     */
    function mulScalar(G1Point memory p, uint256 s) internal view returns (G1Point memory r) {

        uint256[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        // solium-disable-next-line security/no-inline-assembly
        assembly {
            success := staticcall(sub(gas(), 2000), 7, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success, "pairing-mul-failed");
    }
    
    // submod an emelemnt in F_q
    function submod(uint x, uint y) internal pure returns (uint) {
        return addmod(x, PRIME_Q - y, PRIME_Q);
    }

    // Inverse of an element in F_q
    function invmod(uint x) internal pure returns (uint) {
        // Use Euler's theorem
        return expMod(x, PRIME_Q - 2, PRIME_Q);
    }
    // Add two elements in F_q2
    function addmod2(uint[2] memory x, uint[2] memory y) internal pure returns (uint[2] memory) {
        return [addmod(x[0], y[0], PRIME_Q), addmod(x[1], y[1], PRIME_Q)];
    }

    // Subtract two elements in F_q2
    function submod2(uint[2] memory x, uint[2] memory y) internal pure returns (uint[2] memory) {
        return [submod(x[0], y[0]), submod(x[1], y[1])];
    }

    // Multiply two elements in F_q2
    function mulmod2(uint[2] memory x, uint[2] memory y) internal view returns (uint[2] memory) {
        // Use Karatsuba's algorithm
        uint x0y0 = mulmod(x[0], y[0], PRIME_Q);
        uint x1y1 = mulmod(x[1], y[1], PRIME_Q);
        uint x0x1 = addmod(x[0], x[1], PRIME_Q);
        uint y0y1 = addmod(y[0], y[1], PRIME_Q);
        uint x0x1y0y1 = mulmod(x0x1, y0y1, PRIME_Q);
        uint z0 = addmod(x0y0, mulmod(z, x1y1, PRIME_Q), PRIME_Q);
        uint z1 = submod(x0x1y0y1, addmod(x0y0, x1y1, PRIME_Q));
        return [z0, z1];
    }


    // Double a point in G2
    function double(G2Point memory P) internal view returns (G2Point memory) {
        // Use Jacobian coordinates
        if (P.X[0] == 0 && P.X[1] == 0 && P.Y[0] == 0 && P.Y[1] == 0) {
            return P;
        }
        uint[2] memory XX = mulmod2(P.X, P.X);
        uint[2] memory YY = mulmod2(P.Y, P.Y);
        uint[2] memory S = mulmod2([addmod(P.X[0], P.Y[0], PRIME_Q), addmod(P.X[1], P.Y[1], PRIME_Q)], [addmod(P.X[0], P.Y[0], PRIME_Q), addmod(P.X[1], P.Y[1], PRIME_Q)]);
        S = submod2(S, addmod2(XX, YY));
        S = addmod2(S, S);
        uint[2] memory M = addmod2(addmod2(XX, XX), XX);
        uint[2] memory T = submod2(mulmod2(M, M), addmod2(S, S));
        uint[2] memory X3 = T;
        uint[2] memory Y3 = submod2(mulmod2(M, submod2(S, T)), addmod2(YY, YY));
        // uint[2] memory Z3 = addmod2(ZZ, ZZ);
        return G2Point(X3, Y3);
    }

    // Add two points in G2
    function add2(G2Point memory P, G2Point memory Q) internal view returns (G2Point memory) {
        // Use Jacobian coordinates
        if (P.X[0] == 0 && P.X[1] == 0 && P.Y[0] == 0 && P.Y[1] == 0) {
            return Q;
        }
        if (Q.X[0] == 0 && Q.X[1] == 0 && Q.Y[0] == 0 && Q.Y[1] == 0) {
            return P;
        }
        // Use an array to store the intermediate values
        uint[2][12] memory vals;
        vals[0] = mulmod2(P.X, P.X); // Z1Z1
        vals[1] = mulmod2(Q.X, Q.X); // Z2Z2
        vals[2] = mulmod2(P.Y, vals[1]); // U1
        vals[3] = mulmod2(Q.Y, vals[0]); // U2
        vals[4] = mulmod2(mulmod2(P.X, Q.X), vals[1]); // S1
        vals[5] = mulmod2(mulmod2(Q.X, P.X), vals[0]); // S2
        vals[6] = submod2(vals[3], vals[2]); // H
        vals[7] = addmod2(vals[6], vals[6]); // I
        vals[7] = mulmod2(vals[7], vals[7]); // I
        vals[8] = mulmod2(vals[6], vals[7]); // J
        vals[9] = submod2(vals[5], vals[4]); // r
        vals[9] = addmod2(vals[9], vals[9]); // r
        vals[10] = mulmod2(vals[2], vals[7]); // V
        vals[11] = submod2(submod2(mulmod2(vals[9], vals[9]), vals[8]), addmod2(vals[10], vals[10])); // X3
        G2Point memory R; // The result point
        R.X = vals[11];
        R.Y = submod2(mulmod2(vals[9], submod2(vals[10], vals[11])), addmod2(mulmod2(vals[4], vals[8]), vals[4])); // Y3
        return R;
    }

    // A function that performs scalar multiplication on G2 using double-and-add algorithm
    function mulScalar2(G2Point memory p, uint256 s) public view returns (G2Point memory) {
        G2Point memory r = G2Point([uint(0), uint(0)], [uint(0), uint(0)]); // initialize result to zero point
        while (s > 0) {
            if (s % 2 == 1) {
                r = add2(r, p); // add current point if bit is 1
            }
            s = s / 2; // right shift the scalar
            p = double(p); // double the current point
        }
        return r;
    }

    //Return The result of computing the 3-point pairing check
    function pairing_3point(
        G1Point memory a1,
        G2Point memory a2,
        G1Point memory b1,
        G2Point memory b2,
        G1Point memory c1,
        G2Point memory c2
    ) internal view returns (bool) {

        G1Point[3] memory p1 = [a1, b1, c1];
        G2Point[3] memory p2 = [a2, b2, c2];

        uint256 inputSize = 18;
        uint256[] memory input = new uint256[](inputSize);

        for (uint256 i = 0; i < 3; i++) {
            uint256 j = i * 6;
            input[j + 0] = p1[i].X;
            input[j + 1] = p1[i].Y;
            input[j + 2] = p2[i].X[0];
            input[j + 3] = p2[i].X[1];
            input[j + 4] = p2[i].Y[0];
            input[j + 5] = p2[i].Y[1];
        }

        uint256[1] memory out;
        bool success;
        //uint256 len = inputSize * 0x20;
        // solium-disable-next-line security/no-inline-assembly
        assembly {
            success := staticcall(sub(gas(), 2000), 0x8, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success, "pairing-opcode-failed");

        return out[0] != 0;
    }

    /**
     * @dev Returns true if e(a1, a2) * e(b1, b2) == 1
     */
    function pairing_2point(
        G1Point memory a1,
        G2Point memory a2,
        G1Point memory b1,
        G2Point memory b2
    ) internal view returns (bool) {
        // 2 pairs × 6 field elements each = 12
        uint256 inputSize = 12;
        uint256[] memory input = new uint256[](inputSize);

        // first pair
        input[0] = a1.X;
        input[1] = a1.Y;
        input[2] = a2.X[0];
        input[3] = a2.X[1];
        input[4] = a2.Y[0];
        input[5] = a2.Y[1];

        // second pair
        input[6] = b1.X;
        input[7] = b1.Y;
        input[8] = b2.X[0];
        input[9] = b2.X[1];
        input[10] = b2.Y[0];
        input[11] = b2.Y[1];

        uint256[1] memory out;
        bool success;
        // call precompile at address 0x08
        assembly {
            success := staticcall(
                sub(gas(), 2000),
                0x8,
                add(input, 0x20),
                mul(inputSize, 0x20),
                out,
                0x20
            )
            // revert if the call failed
            switch success case 0 { invalid() }
        }
        require(success, "pairing-opcode-failed");
        return out[0] != 0;
    }

    // @dev Modular exponentiation, b^e % _pp.
    function expMod(uint256 _base, uint256 _exp, uint256 _pp) internal pure returns (uint256) {
        require(_pp!=0, "Modulus is zero");

        if (_base == 0)
        return 0;
        if (_exp == 0)
        return 1;

        uint256 r = 1;
        uint256 bit = 57896044618658097711785492504343953926634992332820282019728792003956564819968; // 2 ^ 255
        assembly {
        for { } gt(bit, 0) { }{
            // a touch of loop unrolling for 20% efficiency gain
            r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, bit)))), _pp)
            r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 2))))), _pp)
            r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 4))))), _pp)
            r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 8))))), _pp)
            bit := div(bit, 16)
        }
        }

        return r;
    }

}