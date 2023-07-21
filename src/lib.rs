use std::error::Error;

use ark_bn254::{Bn254, Fr};
use ark_ff::{One, PrimeField};
use ark_groth16::{
    create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
};
use ark_r1cs_std::{
    prelude::{AllocVar, ToBitsGadget},
    R1CSVar,
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError,
};
use ark_serialize::CanonicalSerialize;
use bn::BigUintVar;
use num::{bigint::RandBigInt, BigUint};
use num_prime::RandPrime;
use rand::thread_rng;

mod bn;

const W: usize = 32;
const N: usize = 1024;

struct TestCircuit {
    m: BigUint,
    n: BigUint,
    r: BigUint,
    c: BigUint,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for TestCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let m_var = BigUintVar::<F, W>::new_witness(cs.clone(), || Ok((self.m, N)))?;
        let nn_var = BigUintVar::<F, W>::new_input(cs.clone(), || Ok((&self.n * &self.n, N * 2)))?;
        let g_var =
            BigUintVar::<F, W>::new_input(cs.clone(), || Ok((&self.n + BigUint::one(), N * 2)))?;
        let n_var = BigUintVar::<F, W>::new_input(cs.clone(), || Ok((self.n, N)))?;
        let c_var = BigUintVar::<F, W>::new_input(cs.clone(), || Ok((self.c, N * 2)))?;
        let r_var = BigUintVar::<F, W>::new_witness(cs.clone(), || Ok((self.r, N * 2)))?;

        g_var
            .powm(&m_var.to_bits_le()?, &nn_var, &(BigUint::one() << (N * 2)))?
            .mul_no_carry(&r_var.powm(
                &n_var.to_bits_le()?,
                &nn_var,
                &(BigUint::one() << (N * 2)),
            )?)?
            .rem(&nn_var, &(BigUint::one() << (N * 2)))?
            .enforce_equal_unaligned(&c_var)?;

        Ok(())
    }
}

#[test]
fn test() -> Result<(), Box<dyn Error>> {
    let rng = &mut thread_rng();
    let p: BigUint = rng.gen_prime_exact(N / 2, None);
    let q: BigUint = rng.gen_prime_exact(N / 2, None);
    let n = &p * &q;
    let nn = &n * &n;
    let r = rng.gen_biguint_below(&n);
    let m = rng.gen_biguint_below(&n);
    let g = &n + BigUint::one();

    let cs = ConstraintSystem::new_ref();

    let m_var = BigUintVar::<Fr, W>::new_witness(cs.clone(), || Ok((m.clone(), N)))?;
    let g_var = BigUintVar::<Fr, W>::new_input(cs.clone(), || Ok((g.clone(), N * 2)))?;
    let n_var = BigUintVar::<Fr, W>::new_input(cs.clone(), || Ok((n.clone(), N)))?;
    let nn_var = BigUintVar::<Fr, W>::new_input(cs.clone(), || Ok((nn.clone(), N * 2)))?;
    let r_var = BigUintVar::<Fr, W>::new_witness(cs.clone(), || Ok((r.clone(), N * 2)))?;

    let c_var = g_var
        .powm(&m_var.to_bits_le()?, &nn_var, &(BigUint::one() << (N * 2)))?
        .mul_no_carry(&r_var.powm(&n_var.to_bits_le()?, &nn_var, &(BigUint::one() << (N * 2)))?)?
        .rem(&nn_var, &(BigUint::one() << (N * 2)))?;
    c_var.enforce_lt(&nn_var)?;

    assert_eq!(
        c_var.value()?,
        (g.modpow(&m, &nn) * r.modpow(&n, &nn)) % &nn
    );
    println!("{}", cs.num_constraints());
    assert!(cs.is_satisfied()?);

    Ok(())
}

#[test]
fn test_groth16() -> Result<(), Box<dyn Error>> {
    let rng = &mut thread_rng();
    let p: BigUint = rng.gen_prime_exact(N / 2, None);
    let q: BigUint = rng.gen_prime_exact(N / 2, None);
    let n = &p * &q;
    let nn = &n * &n;
    let r = rng.gen_biguint_below(&n);
    let m = rng.gen_biguint_below(&n);
    let g = &n + BigUint::one();

    let c = (g.modpow(&m, &nn) * r.modpow(&n, &nn)) % &nn;

    let pk = generate_random_parameters::<Bn254, _, _>(
        TestCircuit {
            m: Default::default(),
            n: rng.gen_biguint_range(&(BigUint::one() << N), &(BigUint::one() << N + 1)),
            r: Default::default(),
            c: Default::default(),
        },
        rng,
    )?;
    println!("{}", pk.compressed_size());

    let vk = prepare_verifying_key(&pk.vk);

    let pi = create_random_proof(TestCircuit { m, n: n.clone(), r, c: c.clone() }, &pk, rng)?;

    assert!(verify_proof(
        &vk,
        &pi,
        &vec![
            BigUintVar::<Fr, W>::inputize(&nn, N * 2),
            BigUintVar::<Fr, W>::inputize(&g, N * 2),
            BigUintVar::<Fr, W>::inputize(&n, N),
            BigUintVar::<Fr, W>::inputize(&c, N * 2),
        ]
        .concat()
    )?);

    Ok(())
}
