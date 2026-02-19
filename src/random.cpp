//
//  random.cpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 23/09/2025.
//

#include "random.hpp"

std::mt19937_64 m_mt;
std::uniform_real_distribution<double> rand_uniform;
std::exponential_distribution<double> rand_expo;


double getUni() {
    return rand_uniform(m_mt);
}

double getUniPos() {
    double r;
    do
        r = getUni();
    while ( r * (1. - r) == 0 );
    return r;
}

bool getBool(const double& prob) {
    if (prob == 0.)
        return false;
    else
        return (getUni() <= prob) ? true : false;
}

double getExpo( const double& rate ) {
    
    return rand_expo( m_mt ) / rate ;
    
}

double getErlang( const double& rate, const int& n ) {
    
    double res = 0.;
    for ( int i = 0; i < n; ++i )
        res += getExpo( rate ) ;
    
    return res ;
    
}

/*
 Samples from the equilibrium survivival distribution of X, where X is Erlang.
 Sample k uniformly among 0...n-1. Then draw an Erlang sample with shape n - k
 */
double getErlangSurvival( const double& rate, const int& n ) {

    int k = getUniInt( n - 1 ) ;
    return getErlang( rate, n - k ) ;
    
}

double getGamma(const double& a, const double& b)
{
  /* assume a > 0 */
  int na = floor (a);

  if( a >= std::numeric_limits<int>::max() )
    {
      return b * ( gamma_large ( floor (a) ) + gamma_frac ( a - floor (a) ) );
    }
  else if (a == na)
    {
      return b * gamma_int (na);
    }
  else if (na == 0)
    {
      return b * gamma_frac (a);
    }
  else
    {
      return b * ( gamma_int (na) + gamma_frac (a - na) ) ;
    }
}

double getBeta (const double& a, const double& b)
{
  if ( (a <= 1.0) && (b <= 1.0) )
    {
      double U, V, X, Y;
      while (1)
        {
          U = getUniPos();
          V = getUniPos();
          X = pow( U, 1./a );
          Y = pow( V, 1./b );
          if ( (X + Y) <= 1. )
            {
              if (X + Y > 0)
                {
                  return X/ (X + Y);
                }
              else
                {
                  double logX = log(U)/a;
                  double logY = log(V)/b;
                  double logM = logX > logY ? logX: logY;
                  logX -= logM;
                  logY -= logM;
                  return exp(logX - log(exp(logX) + exp(logY)));
                }
            }
        }
    }
  else
    {
      double x1 = getGamma (a, 1.);
      double x2 = getGamma (b, 1.);
      return x1 / (x1 + x2);
    }
}

int getGeom1( const double& p ) {
    if ( p == 1. )
        return 1;
    return 1 + floor( log( 1. - getUni() ) / log( 1. - p ) );
}

int getBinom(double p, int n) {
    int i, a, b, k = 0;
    while (n > 10) {      /* This parameter is tunable */
        double X;
        a = 1 + (n / 2);
        b = 1 + n - a;
        
        X = getBeta ( (double) a, (double) b );
        
        if (X >= p) {
            n = a - 1;
            p /= X;
        }
        else {
            k += a;
            n = b - 1;
            p = (p - X) / (1 - X);
        }
    }
    
    for (i = 0; i < n; i++) {
        if ( getBool(p) ) {
            k++;
        }
    }
    return k;
}

int getUniInt(const int& n) {
    if (n == 0)
        return 0;
    else
        return static_cast<int>( getUni() * (n + 1) );
} // extrema inclusive


int getPoisson (double mu)
{
  double emu;
  double prod = 1.;
  int k = 0;

  while (mu > 10)
    {
      int m = mu * (7. / 8.);

      double X = gamma_int (m);

      if (X >= mu)
        {
          return k + getBinom (mu / X, m - 1);
        }
      else
        {
          k += m;
          mu -= X;
        }
    }

  /* This following method works well when mu is small */

  emu = exp (-mu);

  do
    {
      prod *= getUni();
      k++;
    }
  while (prod > emu);

  return k - 1;

}

// sample from zero-truncated Poisson distribution
// expected value is mu / ( 1 - e^{-mu} )
int getZeroTruncPoisson( const double &mu ) {
    int res;
    do {
        res = getPoisson( mu );
    }
    while ( res == 0 );
    return res;
}


// sample from negative binomial distribution
// uses same notation as numpy/scipy (NOT wikipedia)
// expected value is mu = n * ( 1 - p ) / p
// variance is var = mu * ( 1 + mu / n )
int getNegBinom( double p, const double n ) {

    if ( p == 1. )
        return 0; //

    double X = getGamma( n, 1. ) ;
    int k = getPoisson( X * ( 1 - p ) / p ) ;
    return k ;
    
}


double gamma_int (const int& a)
{
  if (a < 12)
    {
      int i;
      double prod = 1;

      for (i = 0; i < a; i++)
        {
          prod *= getUniPos();
        }

      /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large ( (double) a );
    }
}


double gamma_large (const double& a)
{
  /* Works only if a > 1, and is most efficient if a is large */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
        {
          y = tan ( M_PI * getUni() );
          x = sqa * y + a - 1;
        }
      while (x <= 0);
      v = getUni();
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

double gamma_frac (const double& a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;

  if (a == 0) {
    return 0;
  }

  p = M_E / (a + M_E);
  do
    {
      u = getUni();
      v = getUniPos();

      if (u < p)
        {
          x = exp ((1 / a) * log (v));
          q = exp (-x);
        }
      else
        {
          x = 1 - log (v);
          q = exp ((a - 1) * log (x));
        }
    }
  while ( getUni() >= q );

  return x;
}
