#pragma once
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;


template <class T>
class TComplex
{
public:
	TComplex();
	TComplex(T real_, T imaginary_);
	TComplex(TComplex<T>& val);
	~TComplex();
	T GetRe();
	T GetIm();
	void SetRe(T val);
	void SetIm(T val);
	T ComplexModule();
	void TrigForm();
	void Degree(double degree);
	void ComplexDegree();
	TComplex operator + (const TComplex<T>& val);
	TComplex operator - (const TComplex<T>& val);
	TComplex operator * (const TComplex<T>& val);
	TComplex operator / (const TComplex<T>& val);
	TComplex& operator += (const TComplex<T>& val);
	TComplex& operator -= (const TComplex<T>& val);
	TComplex& operator *= (const TComplex<T>& val);
	TComplex& operator /= (const TComplex<T>& val);
	TComplex& operator = (const TComplex<T>& val);
	bool operator == (const TComplex<T>& val);
	bool operator != (const TComplex<T>& val);

	template <class D>
	friend ostream& operator << (ostream& out, TComplex<D>& val);
	template <class D>
	friend istream& operator >> (istream& inp, TComplex<D>& val);
	
protected:
	T real, imaginary;
};

template <class T>
inline TComplex<T>::TComplex(): real(1), imaginary(0)
{}

template <class T>
inline TComplex<T>::TComplex(T real_, T imaginary_) : real(real_), imaginary(imaginary_)
{}

template <class T>
inline TComplex<T>::~TComplex()
{}

template <class T>
inline TComplex<T>::TComplex(TComplex& val)
{
	real = val.GetRe();
	imaginary = val.GetIm();
}

template <class T>
inline T TComplex<T>::GetRe()
{
	return real;
}

template <class T>
inline T TComplex<T>::GetIm()
{
	return imaginary;
}

template <class T>
inline void TComplex<T>::SetRe(T val)
{
	real = val;
}

template <class T>
inline void TComplex<T>::SetIm(T val)
{
	imaginary = val;
}

template<class T>
inline T TComplex<T>::ComplexModule()
{
	return (sqrt(real*real + imaginary*imaginary));
}

template<class T>
inline void TComplex<T>::TrigForm()
{
	if ((imaginary == 0) && (real == 0))
	{
		cout << "(" << real << " + " << imaginary << " * i) = "
			<< 0;
	}
	else if (real == 0)
	{
		cout << "(" << real << " + " << imaginary << " * i) = "
			<< imaginary << " * i ";
	}
	else if (imaginary == 0)
	{
		cout << "(" << real << " + " << imaginary << " * i) = "
			<< real;
	}
	else
	{
		T arg = abs((atan(real / imaginary) * 180) / 3.14);
		if ((imaginary > 0) && (real > 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) = "
				<< ComplexModule() << " * " << "cos(" << arg << "+ 2pk) + sin("
				<< arg << "+ 2pk) * i ";
		}
		else if ((imaginary > 0) && (real < 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) = "
				<< ComplexModule() << " * " << "cos(" << -arg << "+ 2pk) + sin("
				<< -arg << "+ 2pk) * i ";
		}
		else if ((imaginary < 0) && (real > 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) = "
				<< ComplexModule() << " * " << "cos(" << 180 - arg << "+ 2pk) + sin("
				<< 180 - arg << "+ 2pk) * i ";
		}
		else
		{
			cout << "(" << real << " + " << imaginary << " * i) = "
				<< ComplexModule() << " * " << "cos(" << arg - 180 << "+ 2pk) + sin("
				<< arg - 180 << "+ 2pk) * i ";
		}
	}
}

template<class T>
inline void TComplex<T>::Degree(double degree)
{
	double IntPart = 0;
	double FracPart = 0;
	FracPart = modf(degree, &IntPart);
	if(FracPart == 0)
	{
		if (degree > 0)
		{
			if ((imaginary == 0) && (real == 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< 0;
			}
			else if (real == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< pow(imaginary, degree) << "i ^ " << degree;
			}
			else if (imaginary == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< pow(real, degree);
			}
			else
			{
				T arg = abs((atan(real / imaginary) * 180) / M_PI);
				if ((imaginary > 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << arg << "+ 2pk) * "
						<< degree << ") + sin((" << arg << "+ 2pk) * " << degree << ") * i )";
				}
				else if ((imaginary > 0) && (real < 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << -arg << "+ 2pk) * "
						<< degree << ") + sin((" << -arg << "+ 2pk) * " << degree << ") * i )";
				}
				else if ((imaginary < 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << 180 - arg << "+ 2pk) * "
						<< degree << ") + sin((" << 180 - arg << "+ 2pk) * " << degree << ") * i )";
				}
				else
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << arg - 180 << "+ 2pk) * "
						<< degree << ") + sin((" << arg - 180 << "+ 2pk) * " << degree << ") * i )";
				}
			}
		}
		else if (degree < 0)
		{
			if ((imaginary == 0) && (real == 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< 0;
			}
			else if (real == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< "(1 / " << pow(ComplexModule(), degree) << " ) * " << "(cos(" << 3.14 / 2 << ") - sin("
					<< 3.14 / 2 << ") * i )";
			}
			else if (imaginary == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< "1 / " << pow(real, degree);
			}
			else
			{
				T arg = abs((atan(real / imaginary) * 180) / 3.14);
				if ((imaginary > 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << arg << "+ 2pk) * "
						<< degree << ") - sin((" << arg << "+ 2pk) * " << degree << ") * i )";
				}
				else if ((imaginary > 0) && (real < 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << -arg << "+ 2pk) * "
						<< degree << ") - sin((" << -arg << "+ 2pk) * " << degree << ") * i )";
				}
				else if ((imaginary < 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << 180 - arg << "+ 2pk) * "
						<< degree << ") - sin((" << 180 - arg << "+ 2pk) * " << degree << ") * i )";
				}
				else
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << arg - 180 << "+ 2pk) * "
						<< degree << ") - sin((" << arg - 180 << "+ 2pk) * " << degree << ") * i )";
				}
			}
		}
		else
		{
			throw("Eror");
		}
	}
	else
	{
		if (degree > 0)
		{
			if ((imaginary == 0) && (real == 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< 0;
			}
			else if (real == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< pow(imaginary, degree) << "i ^ " << degree;
			}
			else if (imaginary == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< pow(real, degree);
			}
			else
			{
				T arg = abs((atan(real / imaginary) * 180) / 3.14);
				if ((imaginary > 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << arg << "+ 2pk) * (1 / "
						<< degree << ")) + sin((" << arg << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
				else if ((imaginary > 0) && (real < 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << -arg << "+ 2pk) * (1 / "
						<< degree << ")) + sin((" << -arg << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
				else if ((imaginary < 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << 180 - arg << "+ 2pk) * "
						<< degree << ") + sin((" << 180 - arg << "+ 2pk) * " << degree << ") * i )";
				}
				else
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
						<< pow(ComplexModule(), degree) << " * " << "(cos((" << arg - 180 << "+ 2pk) * (1 / "
						<< degree << ")) + sin(( " << arg - 180 << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
			}
		}
		else if (degree < 0)
		{
			if ((imaginary == 0) && (real == 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< 0;
			}
			else if (real == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< "(1 / " << pow(ComplexModule(), degree) << " ) * " << "(cos(" << 3.14 / 2 << ") - sin("
					<< 3.14 / 2 << ") * i )";
			}
			else if (imaginary == 0)
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = "
					<< "1 / " << pow(real, degree);
			}
			else
			{
				T arg = abs((atan(real / imaginary) * 180) / 3.14);
				if ((imaginary > 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << arg << "+ 2pk) * (1 / "
						<< degree << ")) - sin((" << arg << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
				else if ((imaginary > 0) && (real < 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << -arg << "+ 2pk) * (1 / "
						<< degree << ")) - sin((" << -arg << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
				else if ((imaginary < 0) && (real > 0))
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << 180 - arg << "+ 2pk) * (1 / "
						<< degree << ")) - sin((" << 180 - arg << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
				else
				{
					cout << "(" << real << " + " << imaginary << " * i) ^ ( " << degree << " ) = 1 / ( "
						<< pow(ComplexModule(), degree) << " ) * " << "(cos((" << arg - 180 << "+ 2pk) * (1 / "
						<< degree << ")) - sin((" << arg - 180 << "+ 2pk) * (1 / " << degree << ")) * i )";
				}
			}
		}
	}
}

template <class T>
inline void TComplex<T>::ComplexDegree()
{
	T c_arg = abs((atan(real / imaginary) * 180) / 3.14);
	TComplex<T> c_degree;
	cin >> c_degree;
	if ((c_degree.real == 0) && (c_degree.imaginary == 0))
	{
		if (((real == 0) && (imaginary == 0)) || ((real == 0) && (imaginary != 0)))
		{
			throw("eror");
		}
		else
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( " << c_degree.real
					 << " + " << c_degree.imaginary << " * i) = " << 0;
		}
	}
	else if(c_degree.real == 0)
	{
		if (((real != 0) && (imaginary != 0)) || ((imaginary != 0) && (real == 0)))
		{
			if ((imaginary > 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< "(1 / (e ^ (" << c_degree.imaginary << " * (" << c_arg << " + 2pk)))"
					<< " * " << "(cos(" << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else if ((imaginary > 0) && (real < 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< "(1 / (e ^ (" << c_degree.imaginary << " * (" << -c_arg << " + 2pk)))"
					<< " * " << "(cos(" << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else if ((imaginary < 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< "(1 / (e ^ (" << c_degree.imaginary << " * (" << 180 - c_arg << " + 2pk)))"
					<< " * " << "(cos(" << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< "(1 / (e ^ (" << c_degree.imaginary << " * (" << c_arg - 180 << " + 2pk)))"
					<< " * " << "(cos(" << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
		}
		else if ((imaginary == 0) && (real != 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< "cos(" << c_degree.imaginary << " * ln(" << real << ") + sin("
				<< c_degree.imaginary << "ln(" << real <<") * i";
		}
		else if ((imaginary == 0) && (real == 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< 1;
		}
	}
	else if (c_degree.imaginary == 0)
	{
		if ((imaginary != 0) && (real != 0))
		{
			if ((imaginary > 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( "
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = 1 / ( "
					<< pow(ComplexModule(), c_degree.real) << " ) * "
					<< "(cos((" << c_arg << "+ 2pk) * (1 / " << c_degree.real << ")) + sin(("
					<< c_arg << "+ 2pk) * (1 / " << c_degree.real << ")) * i )";
			}
			else if ((imaginary > 0) && (real < 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( "
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = 1 / ( "
					<< pow(ComplexModule(), c_degree.real) << " ) * " << "(cos(("
					<< -c_arg << "+ 2pk) * (1 / " << c_degree.real << ")) + sin((" << -c_arg
					<< "+ 2pk) * (1 / " << c_degree.real << ")) * i )";
			}
			else if ((imaginary < 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( "
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = 1 / ( "
					<< pow(ComplexModule(), c_degree.real) << " ) * " << "(cos(("
					<< 180 - c_arg << "+ 2pk) * (1 / " << c_degree.real << ")) + sin((" << 180 - c_arg
					<< "+ 2pk) * (1 / " << c_degree.real << ")) * i )";
			}
			else
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ( "
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = 1 / ( "
					<< pow(ComplexModule(), c_degree.real) << " ) * " << "(cos(("
					<< c_arg - 180 << "+ 2pk) * (1 / " << c_degree.real << ")) + sin((" << c_arg - 180
					<< "+ 2pk) * (1 / " << c_degree.real << ")) * i )";
			}
		}
		else if ((imaginary != 0) && (real == 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< pow(imaginary, c_degree.real) << " * i ^ " << c_degree.real;
		}
		else if ((imaginary == 0) && (real != 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< pow(real, c_degree.real);
		}
		else if ((imaginary == 0) && (real == 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< 1;
		}
	}

	else //(c_degree.real != 0) && (c_degree.imaginary != 0)
	{
		if (((real != 0) && (imaginary != 0)) || ((imaginary != 0) && (real == 0)))
		{
			if ((imaginary > 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< pow(ComplexModule(), c_degree.real) << " * " << "(1 / (e ^ (" << c_degree.imaginary
					<< " * (" << c_arg << " + 2pk)))" << " * " << "(cos(" << c_degree.real << " * ("
					<< c_arg << " + 2pk) + " << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.real << " * (" << c_arg << " + 2pk) + "
					<< c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else if ((imaginary > 0) && (real < 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< pow(ComplexModule(), c_degree.real) << " * " << "(1 / (e ^ (" << c_degree.imaginary
					<< " * (" << -c_arg << " + 2pk)))" << " * " << "(cos(" << c_degree.real << " * ("
					<< -c_arg << " + 2pk) + " << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.real << " * (" << -c_arg << " + 2pk) + "
					<< c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else if ((imaginary < 0) && (real > 0))
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< pow(ComplexModule(), c_degree.real) << " * " << "(1 / (e ^ (" << c_degree.imaginary
					<< " * (" << 180 - c_arg << " + 2pk)))" << " * " << "(cos(" << c_degree.real << " * ("
					<< 180 - c_arg << " + 2pk) + " << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.real << " * (" << 180 - c_arg << " + 2pk) + "
					<< c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
			else
			{
				cout << "(" << real << " + " << imaginary << " * i) ^ ("
					<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
					<< pow(ComplexModule(), c_degree.real) << " * " << "(1 / (e ^ (" << c_degree.imaginary
					<< " * (" << c_arg - 180 << " + 2pk)))" << " * " << "(cos(" << c_degree.real << " * ("
					<< c_arg - 180 << " + 2pk) + " << c_degree.imaginary << " * ln(" << ComplexModule()
					<< ") + (sin(" << c_degree.real << " * (" << c_arg - 180 << " + 2pk) + "
					<< c_degree.imaginary << " * ln(" << ComplexModule() << ") * i";
			}
		}
		else if ((imaginary == 0) && (real != 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< pow(real, c_degree.real) << "cos(" << c_degree.imaginary << " * ln("
				<< real << ") + sin(" << c_degree.imaginary << "ln(" 
				<< real << ") * i";
		}
		else if ((imaginary == 0) && (real == 0))
		{
			cout << "(" << real << " + " << imaginary << " * i) ^ ( "
				<< c_degree.real << " + " << c_degree.imaginary << " * i) = "
				<< 1;
		}
	}
}
template<class T>
inline TComplex<T> TComplex<T>::operator + (const TComplex<T>& val)
{
	return TComplex(real + val.real, imaginary + val.imaginary);
}

template<class T>
inline TComplex<T> TComplex<T>::operator - (const TComplex<T>& val)
{
	return TComplex(real - val.real, imaginary - val.imaginary);
}

template<class T>
inline TComplex<T> TComplex<T>::operator * (const TComplex<T>& val)
{
	return TComplex(real * val.real - imaginary * val.imaginary,
		real * val.imaginary + imaginary * val.real);
}

template<class T>
inline TComplex<T> TComplex<T>::operator / (const TComplex<T>& val)
{
	return TComplex((real * val.real + imaginary * val.imaginary)
		/ (val.real * val.real + val.imaginary * val.imaginary),
		(real * val.imaginary - imaginary * val.real)
		/ (val.real * val.real + val.imaginary * val.imaginary));
}

template<class T>
inline TComplex<T>& TComplex<T>::operator+=(const TComplex<T>& val)
{
	real += val.real;
	imaginary += val.imaginary;
	return*this;
}

template<class T>
inline TComplex<T>& TComplex<T>::operator -= (const TComplex<T>& val)
{
	real -= val.real;
	imaginary -= val.imaginary;
	return *this;
}

template<class T>
inline TComplex<T>& TComplex<T>::operator *= (const TComplex<T>& val)
{
	T re = real;
	T im = imaginary;
	real = re * val.real - im * val.imaginary;
	imaginary = re * val.imaginary + im * val.real;
	return *this;
}

template<class T>
inline TComplex<T>& TComplex<T>::operator /= (const TComplex<T>& val)
{
	T re = real;
	T im = imaginary;
	real = ((re * val.real) + (im * val.imaginary)) / ((val.real * val.real) + (val.imaginary * val.imaginary));
	imaginary = ((im * val.real) - (re * val.imaginary)) / ((val.real * val.real) + (val.imaginary * val.imaginary));
	return *this;
}

template<class T>
inline TComplex<T>& TComplex<T>::operator = (const TComplex<T>& val)
{
	real = val.real;
	imaginary = val.imaginary;
	return*this;
}

template<class T>
inline bool TComplex<T>::operator == (const TComplex<T>& val)
{
	return ((this->real == val.real) && (this->imaginary == val.imaginary));
}

template<class T>
inline bool TComplex<T>::operator != (const TComplex<T>& val)
{
	return ((this->real != val.real) || (this->imaginary != val.imaginary));
}


template <class D>
inline ostream& operator << (ostream& out, TComplex<D>& val)
{
	out << "Complex number = " << val.GetRe() << " + " << val.GetIm() << " * i " << endl;
	system("pause");
	system("cls");
	return out;
}

template <class D>
inline istream& operator >> (istream& inp, TComplex<D>& val)
{
	cout << "Enter real - ";
	inp >> val.real;
	cout << "Enter imaginary - ";
	inp >> val.imaginary;
	system("pause");
	system("cls");
	return inp;
}