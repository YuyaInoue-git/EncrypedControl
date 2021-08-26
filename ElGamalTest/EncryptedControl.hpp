
#ifndef ENCRYPTEDCONTROL
#define ENCRYPTEDCONTROL

#include <cstdint>
#include <ctime>
#include <tuple>
#include <vector>
#include <utility>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/random.hpp>

namespace mp = boost::multiprecision;
namespace rnd = boost::random;

namespace EncryptedControl
{
	class RandomNumber
	{
		public:
			RandomNumber();
			RandomNumber(mp::cpp_int);
			RandomNumber(mp::cpp_int, mp::cpp_int);
			~RandomNumber(void);
			void SetSeed(std::uint64_t);
			void SetRange(mp::cpp_int, mp::cpp_int);
			mp::cpp_int GetRandomNumber(void);

		private:
			mp::cpp_int min;
			mp::cpp_int max;
			std::uint64_t seed;
			rnd::mt19937_64 gen;
			rnd::uniform_int_distribution<mp::cpp_int> dist;
	};

	class ElGamal
	{
		public:
			typedef struct CipherText{
				mp::cpp_int c1;
				mp::cpp_int c2;
			} CipherText;

			typedef struct PublicKey{
				mp::cpp_int p;
				mp::cpp_int q;
				mp::cpp_int g;
				mp::cpp_int h;
			} PublicKey;

			typedef struct SecretKey{
				mp::cpp_int s;
			} SecretKey;

			typedef mp::cpp_dec_float_100 ScalingParameter;

			ElGamal();
			ElGamal(unsigned int);
			~ElGamal();
			std::tuple<PublicKey, SecretKey> KeyGen(unsigned int);
			CipherText Encrypt(mp::cpp_int, PublicKey);
			mp::cpp_int Decrypt(CipherText, SecretKey);
			mp::cpp_int Encode(double, mp::cpp_dec_float_100, PublicKey);
			double Decode(mp::cpp_int, mp::cpp_dec_float_100, PublicKey);
			CipherText Enc(double, mp::cpp_dec_float_100, PublicKey);
			double Dec(CipherText, mp::cpp_dec_float_100, PublicKey, SecretKey);
			std::vector<std::vector<CipherText>> EncMat(std::vector<std::vector<double>>, mp::cpp_dec_float_100, PublicKey);
			std::vector<CipherText> EncVec(std::vector<double>, mp::cpp_dec_float_100, PublicKey);
			std::vector<std::vector<double>> DecMat(std::vector<std::vector<CipherText>>, mp::cpp_dec_float_100, PublicKey, SecretKey);
			std::vector<double> DecVec(std::vector<CipherText>, mp::cpp_dec_float_100, PublicKey, SecretKey);
			std::vector<double> DecMatPlus(std::vector<std::vector<CipherText>>, mp::cpp_dec_float_100, PublicKey, SecretKey);
			PublicKey GetPublicKey(void);
			SecretKey GetSecretKey(void);

		private:
			RandomNumber randomNumber;
			PublicKey pk;
			SecretKey sk;
			unsigned int keyLength;

			bool IsElement(mp::cpp_int, mp::cpp_int, mp::cpp_int);
			bool IsGenerator(mp::cpp_int, mp::cpp_int, mp::cpp_int);
			mp::cpp_int GetGenerator(mp::cpp_int, mp::cpp_int);
	};

	class ModUtil
	{
		public:
			static mp::cpp_int Mod(mp::cpp_int a, mp::cpp_int m)
			{
				mp::cpp_int c;

				c = a % m;

				return c < 0? c + mp::abs(m) : c;
			}

			static mp::cpp_int ModPow(mp::cpp_int a, mp::cpp_int b, mp::cpp_int m)
			{
				mp::cpp_int c;

				c = mp::powm(a, b, m);

				return c < 0? c + mp::abs(m) : c;
			}

			static mp::cpp_int ModMul(mp::cpp_int a, mp::cpp_int b, mp::cpp_int m)
			{
				mp::cpp_int c;

				c = (a * b) % m;

				return c < 0? c + mp::abs(m) : c;
			}

			static mp::cpp_int ModInv(mp::cpp_int a, mp::cpp_int m)
			{
				mp::cpp_int b, u, v, t;

				b = m;
				u = 1;
				v = 0;

				while(b != 0){
					t = a / b;
					a -= t * b;
					u -= t * v;
					std::swap(a, b);
					std::swap(u, v);
				}

				u %= m;
				return u < 0? u + mp::abs(m) : u;
			}

			static ElGamal::CipherText ModHadamardProduct(ElGamal::CipherText A, ElGamal::CipherText B, mp::cpp_int m)
			{
				ElGamal::CipherText C;

				C.c1 = ModMul(A.c1, B.c1, m);
				C.c2 = ModMul(A.c2, B.c2, m);

				return C;
			}
	};

	class PrimeUtil
	{
		public:
			static bool IsPrime(mp::cpp_int q)
			{
				rnd::mt19937_64 gen;

				return mp::miller_rabin_test(q, 50, gen);
			}

			static mp::cpp_int GetSafePrime(unsigned int bitLength)
			{
				mp::cpp_int binMin = static_cast<mp::cpp_int>(1) << (bitLength - 1);
				mp::cpp_int binMax = binMin << 1;

				RandomNumber randomNumber(binMin, binMax - 1);

				mp::cpp_int q = randomNumber.GetRandomNumber();
				q = mp::bit_set(q, 0);
				q = mp::bit_set(q, bitLength - 1);

				while(PrimeUtil::IsPrime(q) == false || IsPrime(2 * q + 1) == false){
					q = randomNumber.GetRandomNumber();
					q = mp::bit_set(q, 0);
					q = mp::bit_set(q, bitLength - 1);
				}

				return q;
			}
	};

	class EncryptedController
	{
		public:
			static std::vector<std::vector<ElGamal::CipherText>> GetSignal(std::vector<std::vector<ElGamal::CipherText>> C_Phi, std::vector<ElGamal::CipherText> C_xi, ElGamal::PublicKey pk)
			{
				std::vector<std::vector<ElGamal::CipherText>> C_Psi(C_Phi.size(), std::vector<ElGamal::CipherText>(C_Phi[0].size()));

				for(int i = 0; i < C_Phi.size(); i ++){
					for(int j = 0; j < C_Phi[0].size(); j ++){
						C_Psi[i][j] = ModUtil::ModHadamardProduct(C_Phi[i][j], C_xi[j], pk.p);
					}
				}

				return C_Psi;
			}
	};
}

using namespace EncryptedControl;

RandomNumber::RandomNumber()
	: min(0), max(1)
{
	RandomNumber::SetSeed(static_cast<std::uint64_t>(std::time(0)));
	RandomNumber::SetRange(min, max);
}

RandomNumber::RandomNumber(mp::cpp_int a)
	: min(0), max(a)
{
	RandomNumber::SetSeed(static_cast<std::uint64_t>(std::time(0)));
	RandomNumber::SetRange(min, max);
}

RandomNumber::RandomNumber(mp::cpp_int a, mp::cpp_int b)
	: min(a), max(b)
{
	RandomNumber::SetSeed(static_cast<std::uint64_t>(std::time(0)));
	RandomNumber::SetRange(min, max);
}

RandomNumber::~RandomNumber()
{
	
}

void RandomNumber::SetSeed(std::uint64_t new_seed)
{
	gen.seed(new_seed);
}

void RandomNumber::SetRange(mp::cpp_int min, mp::cpp_int max)
{
	rnd::uniform_int_distribution<mp::cpp_int>::param_type newParam(min, max);
	dist.param(newParam);
}

mp::cpp_int RandomNumber::GetRandomNumber()
{
	return dist(gen);
}

ElGamal::ElGamal()
	: keyLength(10)
{
	std::tie(pk, sk) = KeyGen(keyLength);
}

ElGamal::ElGamal(unsigned int k)
	: keyLength(k)
{
	std::tie(pk, sk) = KeyGen(keyLength);
}

ElGamal::~ElGamal()
{

}

std::tuple<ElGamal::PublicKey, ElGamal::SecretKey> ElGamal::KeyGen(unsigned int keyLength)
{
	ElGamal::PublicKey pk;
	ElGamal::SecretKey sk;

	pk.q = PrimeUtil::GetSafePrime(keyLength);
	pk.p = 2 * pk.q + 1;

	pk.g = ElGamal::GetGenerator(pk.q, pk.p);

	randomNumber.SetRange(0, pk.q - 1);

	sk.s = randomNumber.GetRandomNumber();
	while(sk.s == 0){
		sk.s = randomNumber.GetRandomNumber();
	}

	pk.h = ModUtil::ModPow(pk.g, sk.s, pk.p);

	return std::forward_as_tuple(pk, sk);
}

ElGamal::CipherText ElGamal::Encrypt(mp::cpp_int m, ElGamal::PublicKey pk)
{
	ElGamal::CipherText C;

	mp::cpp_int r = randomNumber.GetRandomNumber();

	C.c1 = ModUtil::ModPow(pk.g, r, pk.p);
	C.c2 = ModUtil::ModMul(m, ModUtil::ModPow(pk.h, r, pk.p), pk.p);

	return C;
}

mp::cpp_int ElGamal::Decrypt(ElGamal::CipherText C, ElGamal::SecretKey sk)
{
	mp::cpp_int m;

	m = ModUtil::ModMul(C.c2, ModUtil::ModInv(ModUtil::ModPow(C.c1, sk.s, pk.p), pk.p), pk.p);

	return m;
}

mp::cpp_int ElGamal::Encode(double x, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk)
{
	mp::cpp_int m;
	mp::cpp_int firstDecimalPlace;
	mp::cpp_int upperCandidate;
	mp::cpp_int lowerCandidate;

	m = static_cast<mp::cpp_int>(static_cast<mp::cpp_dec_float_100>(x) * gamma);
	firstDecimalPlace = ModUtil::Mod(static_cast<mp::cpp_int>(static_cast<mp::cpp_dec_float_100>(x) * gamma * 10), 10);

	if(m < 0){
		if(m < -(pk.q + 1)){
			return -1;
		}
		else{
			m += pk.p;
		}
	}
	else{
		if(m >= pk.q){
			return -1;
		}
	}

	if(ElGamal::IsElement(m, pk.q, pk.p) == true){
		return m;
	}
	else{
		lowerCandidate = m - 1;
		upperCandidate = m + 1;
		
		if(firstDecimalPlace >= 5 || firstDecimalPlace == 0){
			for(int i = 0; i < pk.q; i ++){
				if(ElGamal::IsElement(upperCandidate, pk.q, pk.p) == true){
					return upperCandidate;
				}
				else if(ElGamal::IsElement(lowerCandidate, pk.q, pk.p) == true){
					return lowerCandidate;
				}
				else{
					lowerCandidate --;
					upperCandidate ++;
				}
			}
		}
		else{
			for(int i = 0; i < pk.q; i ++){
				if(ElGamal::IsElement(lowerCandidate, pk.q, pk.p) == true){
					return lowerCandidate;
				}
				else if(ElGamal::IsElement(upperCandidate, pk.q, pk.p) == true){
					return upperCandidate;
				}
				else{
					lowerCandidate --;
					upperCandidate ++;
				}
			}
		}
	}

	return -1;
}

double ElGamal::Decode(mp::cpp_int m, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk)
{
	if(m > pk.q){
		m -= pk.p;
	}

	return static_cast<double>(static_cast<mp::cpp_dec_float_100>(m) / gamma);
}

ElGamal::CipherText ElGamal::Enc(double x, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk)
{
	return ElGamal::Encrypt(ElGamal::Encode(x, gamma, pk), pk);
}

double ElGamal::Dec(ElGamal::CipherText C, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk, ElGamal::SecretKey sk)
{
	return ElGamal::Decode(ElGamal::Decrypt(C, sk), gamma, pk);
}

std::vector<std::vector<ElGamal::CipherText>> ElGamal::EncMat(std::vector<std::vector<double>> A, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk)
{
	std::vector<std::vector<ElGamal::CipherText>> C(A.size(), std::vector<ElGamal::CipherText>(A[0].size()));

	for(int i = 0; i < A.size(); i ++){
		for(int j = 0; j < A[0].size(); j ++){
			C[i][j] = ElGamal::Enc(A[i][j], gamma, pk);
		}
	}

	return C;
}

std::vector<ElGamal::CipherText> ElGamal::EncVec(std::vector<double> v, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk)
{
	std::vector<ElGamal::CipherText> C(v.size());

	for(int i = 0; i < v.size(); i ++){
		C[i] = ElGamal::Enc(v[i], gamma, pk);
	}

	return C;
}

std::vector<std::vector<double>> ElGamal::DecMat(std::vector<std::vector<ElGamal::CipherText>> C, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk, ElGamal::SecretKey sk)
{
	std::vector<std::vector<double>> A(C.size(), std::vector<double>(C[0].size()));

	for(int i = 0; i < C.size(); i ++){
		for(int j = 0; j < C[0].size(); j ++){
			A[i][j] = ElGamal::Dec(C[i][j], gamma, pk, sk);
		}
	}

	return A;
}

std::vector<double> ElGamal::DecVec(std::vector<ElGamal::CipherText> C, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk, ElGamal::SecretKey sk)
{
	std::vector<double> v(C.size());

	for(int i = 0; i < C.size(); i ++){
		v[i] = ElGamal::Dec(C[i], gamma, pk, sk);
	}

	return v;
}

std::vector<double> ElGamal::DecMatPlus(std::vector<std::vector<ElGamal::CipherText>> C, mp::cpp_dec_float_100 gamma, ElGamal::PublicKey pk, ElGamal::SecretKey sk)
{
	std::vector<double> v(C.size(), 0);

	for(int i = 0; i < C.size(); i ++){
		for(int j = 0; j < C[0].size(); j ++){
			v[i] += ElGamal::Dec(C[i][j], gamma, pk, sk);
		}
	}

	return v;
}

bool ElGamal::IsElement(mp::cpp_int m, mp::cpp_int q, mp::cpp_int p)
{
	if(ModUtil::ModPow(m, q, p) == 1){
		return true;
	}
	else{
		return false;
	}
}

bool ElGamal::IsGenerator(mp::cpp_int g, mp::cpp_int q, mp::cpp_int p)
{
	if(ModUtil::ModPow(g, q, p) == 1){
		return true;
	}
	else{
		return false;
	}
}

mp::cpp_int ElGamal::GetGenerator(mp::cpp_int q, mp::cpp_int p)
{
	mp::cpp_int g = 2;

	while(IsGenerator(g, q, p) == false){
		g ++;
	}

	if(g < p){
		return g;
	}
	else{
		return -1;
	}
}

ElGamal::PublicKey ElGamal::GetPublicKey()
{
	return pk;
}

ElGamal::SecretKey ElGamal::GetSecretKey()
{
	return sk;
}

#endif