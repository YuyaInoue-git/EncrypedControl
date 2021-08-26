// ElGamalTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "EncryptedControl.hpp"

using namespace EncryptedControl;

int main()
{
	/* real number */
	double x1, x2, x3;

	/* plaintext */
	mp::cpp_int m1, m2, m3;

	/* ciphertext */
	ElGamal::CipherText C1, C2, C3;

	/* key length (bit) */
	unsigned int keyLength = 64;

	/* scaling parameter */
	ElGamal::ScalingParameter gamma = 1e8;

	ElGamal cryptosystem(keyLength);

	x1 = 1.23456789;
	x2 = 0.0123456789;

	/* encode */
	m1 = cryptosystem.Encode(x1, gamma, cryptosystem.GetPublicKey());
	m2 = cryptosystem.Encode(x2, gamma, cryptosystem.GetPublicKey());

	/* encrypt */
	C1 = cryptosystem.Encrypt(m1, cryptosystem.GetPublicKey());
	C2 = cryptosystem.Encrypt(m2, cryptosystem.GetPublicKey());

	/* Hadamard product */
	C3 = ModUtil::ModHadamardProduct(C1, C2, cryptosystem.GetPublicKey().p);

	/* decrypt */
	m3 = cryptosystem.Decrypt(C3, cryptosystem.GetSecretKey());

	/* decode */
	x3 = cryptosystem.Decode(m3, gamma * gamma, cryptosystem.GetPublicKey());

	std::cout << "true = " << std::fixed << std::setprecision(10) << x1 * x2 << std::endl;
	std::cout << "result = " << std::fixed << std::setprecision(10) << x3 << std::endl;

	return 0;
}
