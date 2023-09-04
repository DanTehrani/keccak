#![feature(int_roundings)]

// Keccak256 parameters in bits
const ROUNDS: usize = 24;
const OUTPUT_LEN: usize = 256;
const CAPACITY: usize = OUTPUT_LEN * 2;
const STATE_WIDTH: usize = 1600;
const RATE: usize = STATE_WIDTH - CAPACITY;

// Table 2 of https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.202.pdf
const RHO_OFFSETS: [[u32; 5]; 5] = [
    [0, 1, 190, 28, 91],
    [36, 300, 6, 55, 276],
    [3, 10, 171, 153, 231],
    [105, 45, 15, 21, 136],
    [210, 66, 253, 120, 78],
];

// Copied from https://github.com/debris/tiny-keccak/blob/master/src/keccakf.rs
const RC: [u64; ROUNDS] = [
    1u64,
    0x8082u64,
    0x800000000000808au64,
    0x8000000080008000u64,
    0x808bu64,
    0x80000001u64,
    0x8000000080008081u64,
    0x8000000000008009u64,
    0x8au64,
    0x88u64,
    0x80008009u64,
    0x8000000au64,
    0x8000808bu64,
    0x800000000000008bu64,
    0x8000000000008089u64,
    0x8000000000008003u64,
    0x8000000000008002u64,
    0x8000000000000080u64,
    0x800au64,
    0x800000008000000au64,
    0x8000000080008081u64,
    0x8000000000008080u64,
    0x80000001u64,
    0x8000000080008008u64,
];

fn to_words(bytes: &[u8]) -> Vec<u64> {
    let mut bytes_resized = bytes.to_vec();
    bytes_resized.resize(bytes.len().next_multiple_of(8), 0);

    let mut words = vec![0u64; bytes_resized.len() / 8];
    for i in 0..words.len() {
        words[i] = u64::from_le_bytes([
            bytes_resized[i * 8],
            bytes_resized[i * 8 + 1],
            bytes_resized[i * 8 + 2],
            bytes_resized[i * 8 + 3],
            bytes_resized[i * 8 + 4],
            bytes_resized[i * 8 + 5],
            bytes_resized[i * 8 + 6],
            bytes_resized[i * 8 + 7],
        ]);
    }
    words
}

// Multi-rate padding as described in Section 5 of https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.202.pdf
fn pad(n_bytes: usize) -> Vec<u8> {
    let padded_size_in_bytes = n_bytes.next_multiple_of(RATE / 8);
    let mut padding = vec![0u8; padded_size_in_bytes - n_bytes];
    let padding_len = padding.len();
    padding[0] = 0x01;
    padding[padding_len - 1] = 0x80;
    padding
}

pub fn keccak256(input: &[u8]) -> [u8; 32] {
    if input.len() > (RATE / 8) {
        panic!("Sponge not supported");
    };

    // Padding
    let pad_vec = pad(input.len());
    let mut p = input.to_vec();
    p.extend(pad_vec);

    // Transform padded input into array of words
    let mut words = to_words(&p);
    words.resize(25, 0);
    let mut state: [u64; 25] = words.try_into().unwrap();

    // The permutation
    keccak_p(&mut state);

    // Output in bytes (32 bytes)
    let out = [state[0], state[1], state[2], state[3]];
    let mut out_bytes = [0u8; 32];
    out_bytes[0..8].copy_from_slice(&out[0].to_le_bytes());
    out_bytes[8..16].copy_from_slice(&out[1].to_le_bytes());
    out_bytes[16..24].copy_from_slice(&out[2].to_le_bytes());
    out_bytes[24..32].copy_from_slice(&out[3].to_le_bytes());

    out_bytes
}

pub fn keccak_p(state: &mut [u64; 25]) {
    for i in 0..ROUNDS {
        // ############################################
        // Theta
        // ############################################

        let mut c = [0; 5];
        let mut d = [0; 5];

        for y in 0..5 {
            for x in 0..5 {
                c[x] ^= state[x + y * 5];
            }
        }

        for x in 0..5 {
            d[x] = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
        }

        for y in 0..5 {
            for x in 0..5 {
                state[x + y * 5] ^= d[x];
            }
        }

        // ############################################
        // Rho
        // ############################################

        let mut rho_x = 0;
        let mut rho_y = 1;
        for _ in 0..24 {
            // Rotate each lane by an offset
            let index = rho_x + 5 * rho_y;
            state[index] = state[index].rotate_left(RHO_OFFSETS[rho_y][rho_x] % 64);
            let rho_x_prev = rho_x;
            rho_x = rho_y;
            rho_y = (2 * rho_x_prev + 3 * rho_y) % 5;
        }

        // ############################################
        // Pi
        // ############################################

        let state_cloned = state.clone();
        for y in 0..5 {
            for x in 0..5 {
                let index = ((x + 3 * y) % 5) + x * 5;
                state[x + y * 5] = state_cloned[index];
            }
        }

        // ############################################
        // Chi
        // ############################################

        let state_cloned = state.clone();
        for y in 0..5 {
            for x in 0..5 {
                let index = x + y * 5;
                state[index] = state_cloned[index]
                    ^ (!state_cloned[(x + 1) % 5 + y * 5]) & state_cloned[(x + 2) % 5 + y * 5];
            }
        }

        // ############################################
        // Iota
        // ############################################

        state[0] ^= RC[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use tiny_keccak::{Hasher, Keccak};

    #[test]
    fn test_keccak_p() {
        let input = (0..32).collect::<Vec<u8>>();

        // Compare with tiny-keccak
        let mut ref_hasher = Keccak::v256();
        ref_hasher.update(&input);
        let mut expected = [0u8; 32];
        ref_hasher.finalize(&mut expected);

        let result = keccak256(&input);
        assert_eq!(&result, &expected);
    }
}
