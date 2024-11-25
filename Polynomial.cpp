#include "Polynomial.h"

Polynomial::Polynomial() : coeffs(1, 0.0) {}

Polynomial::Polynomial(const vector<double>& coefficients) : coeffs(coefficients) {
  // Remove trailing zeros
  while (!coeffs.empty() && coeffs.back() == 0.0) {
    coeffs.pop_back();
  }
}

Polynomial::Polynomial(const Polynomial& other) : coeffs(other.coeffs) {}

Polynomial::~Polynomial() {}

Polynomial& Polynomial::operator=(const Polynomial& other) {
  coeffs = other.coeffs;
  return *this;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
  int maxDegree = std::max(degree(), other.degree());
  std::vector<double> result(maxDegree + 1, 0.0);
  for (int i = 0; i <= degree(); ++i) {
    result[i] += coeffs[i];
  }
  for (int i = 0; i <= other.degree(); ++i) {
    result[i] += other.coeffs[i];
  }
  return Polynomial(result);
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
  int maxDegree = std::max(degree(), other.degree());
  std::vector<double> result(maxDegree + 1, 0.0);
  for (int i = 0; i <= degree(); ++i) {
    result[i] += coeffs[i];
  }
  for (int i = 0; i <= other.degree(); ++i) {
    result[i] -= other.coeffs[i];
  }
  return Polynomial(result);
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
  int n = degree() + other.degree();
  std::vector<double> result(n + 1, 0.0);
  for (int i = 0; i <= degree(); ++i) {
    for (int j = 0; j <= other.degree(); ++j) {
      result[i + j] += coeffs[i] * other.coeffs[j];
    }
  }
  return Polynomial(result);
}

bool Polynomial::operator==(const Polynomial& other) const {
  return coeffs == other.coeffs;
}

std::ostream& operator<<(std::ostream& out, const Polynomial& poly) {
  if (poly.coeffs.empty()) {
    out << "0";
    return out;
  }
  bool firstTerm = true;
  for (int i = poly.degree(); i >= 0; --i) {
    if (poly.coeffs[i] != 0.0) {
      if (!firstTerm) {
        out << (poly.coeffs[i] > 0 ? " + " : " - ");
      }
      firstTerm = false;
      if (std::abs(poly.coeffs[i]) != 1 || i != 0) {
        out << std::abs(poly.coeffs[i]);
      }
      if (i > 0) {
        out << "x^" << i;
      }
    }
  }
  return out;
}

int Polynomial::degree() const {
  return coeffs.size() - 1;
}

double Polynomial::evaluate(double x) const {
  double result = 0.0;
  for (int i = degree(); i >= 0; --i) {
    result += coeffs[i] * std::pow(x, i);
  }
  return result;
}

Polynomial Polynomial::derivative() const {
  if (coeffs.empty()) {
    return Polynomial();
  }
  std::vector<double> result(degree());
  for (int i = 1; i <= degree(); ++i) {
    result[i - 1] = coeffs[i] * i;
  }
  return Polynomial(result);
}

Polynomial Polynomial::integral() const {
  if (coeffs.empty()) {
    return Polynomial();
  }
  std::vector<double> result(degree() + 2, 0.0);
  for (int i = 0; i <= degree(); ++i) {
    result[i + 1] = coeffs[i] / (i + 1);
  }
  return Polynomial(result);
}

double Polynomial::integral(double x1, double x2) const {
  Polynomial integralPoly = integral();
  return integralPoly.evaluate(x2) - integralPoly.evaluate(x1);
}

double Polynomial::getRoot(double guess, double tolerance, int maxIter) const {
  for (int i = 0; i < maxIter; ++i) {
    double f = evaluate(guess);
    double df = derivative().evaluate(guess);
    if (std::abs(df) < 1e-10) {
      break; // Avoid division by zero
    }
    double newGuess = guess - f / df;
    if (std::abs(newGuess - guess) < tolerance) {
      return newGuess;
    }
    guess = newGuess;
  }
  return guess; // Return the last guess if convergence is not achieved
}

void Polynomial::setCoefficients(const vector<double>& coefficients) {
  coeffs = coefficients;
}

double Polynomial::getCoefficient(int degree) const {
  if (degree < 0 || degree > this->degree()) {
    return 0.0;
  }
  return coeffs[degree];
}