/*************************************************/
/*************************************************/
/* */
/* Decomposition of FIR Filter */
/* */
/*************************************************/
/*************************************************/
#include <complex.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NPLOT 2048
#define PI 3.14159265358979363846
class Polynomial
{
public:
  void read_from_file();
  unsigned char get_roots();
  void write_to_file();
  void write_to_screen();
  void write_error_message(unsigned char);
  void check_value();
  ~Polynomial();

private:
  unsigned char poly_check();
  void quadratic(std::complex<double> *);
  unsigned char lin_or_quad(std::complex<double> *);
  void hornc(std::complex<double>, unsigned char);
  void horncd(double, double);
  int poldef(unsigned char);
  void monic();
  // functions blow are for Newton's Method
  std::complex<double> newton(std::complex<double>, double *);
  void f_value1(std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double>);
  void f_value2(std::complex<double> *, std::complex<double> *, std::complex<double>);
  // functions blow are for Muller's Method

  std::complex<double> muller();
  void initialize(std::complex<double> *, double *);
  void root_of_parabola();
  void iteration_equation(double *);
  void suppress_overflow();
  void too_big_functionvalues(double *);
  void convergence_check(int *, double, double, double);
  void compute_function(double, double *, double);
  void check_x_value(std::complex<double> *, double *, int *, double, double, double, int *);
  void root_check(double, int *, int *, int *, std::complex<double>);
  void f_value(int, std::complex<double> *, std::complex<double>);
  std::complex<double> x0, x1, x2, // common points [x0,f(x0)=P(x0)], ... [x2,f(x2)]
      f0, f1, f2,                  // of parabola and polynomial
      h1, h2,                      // distance between x2 and x1
      q2, *psave, *psave1;         // smaller root of parabola
  int iter,
      nred, // the highest exponent of the deflated polynomial
      n, N; // original degree of the input
  int distinct, indicator;
  double data;
  double maxerr;
  double *matlab;
  static std::complex<double> *p, // coefficient vector of polynomial
      *pred,                      // coefficient vector of deflated polynom.
      *root;                      // vector of determined roots
  static unsigned char flag;
};
unsigned char Polynomial::flag = 1;
std::complex<double> *Polynomial::p;
std::complex<double> *Polynomial::pred;
std::complex<double> *Polynomial::root;
// read coefficients stored in file FILENAME
void Polynomial ::read_from_file()
{
  char filename[32];
  int i; // counter
  FILE *file_ptr;
  // open file

  printf("*******************************\n");
  printf("*** FIR Filter Decompostion ***\n");
  printf("*******************************\n\n");
  printf("Enter filter file name: ");
  scanf("%s", filename);
  if ((file_ptr = fopen(filename, "r")) == NULL)
  {
    printf("Can't open file %s!\n", filename);
    exit(0);
  }
  else
  {
    file_ptr = fopen(filename, "r");
    // read degree
    (void)fscanf(file_ptr, "%d %d %d", &n, &distinct, &indicator);
    // allocate the memory
    N = n;
    p = new std::complex<double>[n];
    pred = new std::complex<double>[n];
    root = new std::complex<double>[n-1];
    psave1 = p;
    psave = pred;
    flag = 0;
    // read coefficients
    for (int i = n - 1; i >= distinct - 1; i--)
    {
      fscanf(file_ptr, "%lf", &data);
      if (fabs(data) < 1e-8)
        data = 0.0;
      p[i].real(data);
      p[i].imag(0.0);
    }
    if (indicator == 0)
    {
      for (i = 0; i < distinct; i++)
      {
        p[i].real(p[n - i - 1].real());
        p[i].imag(0.0);
      }
    }
    else
    {
      for (i = 0; i < distinct; i++)
      {
        p[i].real(-p[n - i - 1].real());
        p[i].imag(0.0);
      }
      if (n % 2 != 0)
        p[distinct - 1].real(-p[distinct - 1].real());
    }
    (void)fclose(file_ptr);
    file_ptr = fopen("coeffs_out.dat", "w");
    fprintf(file_ptr, "[\n");

    for (i = 0; i < n; i++)
    {
      fprintf(file_ptr, "%.8f,\n", p[i].real());
    }


    fprintf(file_ptr, "]\n");
    (void)fclose(file_ptr);

  }
  matlab = new double[N];
  for (i = n - 1; i >= 0; i--)
    matlab[i] = p[i].real();
}
// write the roots result to a file
void Polynomial ::write_to_file()
{
  char sig0, sig1,      // sign of real part
      sig2, sig3, sig4; // sign of imaginary part
  int i, k = 0;
  int count = 0;
  int out_uc = 0;
  int in_uc = 0;
  int on_uc = 0;
  int on_axis = 0;
  int on_one = 0;
  double *h_0, *h_1, *h_2, *h_3, *h_4;
  double *temp_0, *temp_1, *temp_2, *temp_3, *temp_4;
  FILE *file_ptr;
  std::complex<double> *inside;
  std::complex<double> *outside;
  std::complex<double> *onunitcircle;
  std::complex<double> *onaxis;
  std::complex<double> *onone;
  std::complex<double> Hm, temp;
  Hm.real(1.0);
  Hm.imag(1.0);
  // allocate the memory
  inside = new std::complex<double>[n];
  outside = new std::complex<double>[n];
  onunitcircle = new std::complex<double>[n];
  onaxis = new std::complex<double>[n];
  onone = new std::complex<double>[n];
  h_0 = new double[n];
  h_1 = new double[n];
  h_2 = new double[n];
  h_3 = new double[n];
  h_4 = new double[n];
  temp_0 = new double[n];
  temp_1 = new double[n];
  temp_2 = new double[n];
  temp_3 = new double[n];
  temp_4 = new double[n];
  for (i = 0; i < n; i++)
  {
    h_0[i] = 1.0;
    h_1[i] = 1.0;
    h_2[i] = 1.0;
    h_3[i] = 0.0;
    h_4[i] = 0.0;
  }
  // generate output file
  file_ptr = fopen("roots_FIR.dat", "w");
  fprintf(file_ptr, " %6d \n", n);
  for (i = 0; i < n; i++)
  {
    sig1 = (root[i].real() >= 0) ? ' ' : '-';
    sig2 = (root[i].imag() >= 0) ? ' ' : '-';
    fprintf(file_ptr, " %c%.18e %c%.18e\n", sig1, fabs(root[i].real()), sig2,
            fabs(root[i].imag()));
  }
  (void)fclose(file_ptr);
  file_ptr = fopen("roots_FIR.txt", "w");
  fprintf(file_ptr, " Decomposition of FIR Filter\n");
  fprintf(file_ptr, " Roots of the FIR Filter\n");
  fprintf(file_ptr, " Filter Order = %d\n\n", n);
  fprintf(file_ptr, " Real Part Imaginary Part\n");
  for (i = 0; i < n; i++)
  {
    sig1 = (root[i].real() >= 0) ? ' ' : '-';
    sig2 = (root[i].imag() >= 0) ? ' ' : '-';
    fprintf(file_ptr, " %c%.18e %c%.18e\n", sig1, fabs(root[i].real()), sig2,
            fabs(root[i].imag()));
  }
  (void)fclose(file_ptr);
  
  
  for (i = 0; i < n; i++)
  {
    if (sqrt(root[i].real() * root[i].real() + root[i].imag() * root[i].imag()) < (1 - 1E-8) &&
        (root[i].imag() > 0))
    {
      inside[in_uc].real(root[i].real());
      inside[in_uc].imag(root[i].imag());
      in_uc++;
    }
  }
  // save zeroes outside the unit circle. (Above - Half)
  for (i = 0; i < n; i++)
  {
    if (sqrt(root[i].real() * root[i].real() + root[i].imag() * root[i].imag()) - 1 > (1E-8) &&
        (root[i].imag() > 0))
    {
      outside[out_uc].real(root[i].real());
      outside[out_uc].imag(root[i].imag());
      out_uc++;
    }
  }
  // save zeroes on the unit circle.
  for (i = 0; i < n; i++)
  {
    if (sqrt(root[i].real() * root[i].real() + root[i].imag() * root[i].imag()) < (1 + 1E-8) &&
        sqrt(root[i].real() * root[i].real() + root[i].imag() * root[i].imag()) > (1 - 1E-8) &&
        (fabs(root[i].imag()) != 0))
    {
      onunitcircle[on_uc].real(root[i].real());
      onunitcircle[on_uc].imag(root[i].imag());
      on_uc++;
    }
  }
  // save zeroes on the Real Axis.(All )
  for (i = 0; i < n; i++)
  {
    if (fabs(root[i].imag()) == 0 &&
        ((fabs(root[i].real()) > 1 + 1E-8) || (fabs(root[i].real()) < 1 - 1E-8)))
    {
      onaxis[on_axis].real(root[i].real());
      onaxis[on_axis].imag(root[i].imag());
      on_axis++;
    }
  }
  // save zeroes on the Real Axis are +1 or -1.
  for (i = 0; i < n; i++)
  {
    if (((fabs(root[i].imag()) == 0) && (fabs(root[i].real()) < 1 + 1E-8) &&
         (fabs(root[i].real()) > 1 - 1E-8)))
    {
      onone[on_one].real(root[i].real());
      onone[on_one].imag(root[i].imag());
      on_one++;
    }
  }
  // get the subfilter frequency response
  file_ptr = fopen("FIR_subfilter_response.dat", "w");
  for (i = 0; i < on_one; i++)
  {
    if ((onone[i].real() - 1) > -1E-8 && (onone[i].real() - 1) < 1E-8)
    {
      h_1[i] = -1;
    }
    temp_0[count] = h_0[i];
    temp_1[count] = h_1[i];
    temp_2[count] = temp_3[count] = temp_4[count] = 0.0;
    count++;
  }
  for (i = 0; i < on_axis; i++)
  {
    if (fabs(onaxis[i].real()) > 1)
      k++;
    else if (fabs(onaxis[i].real()) < 1)
    {
      if (onaxis[i].real() != 0 || onaxis[i].imag() != 0)
      {
        h_1[i + on_one - k] = -(onaxis[i].real() + 1 / onaxis[i].real());
        temp_0[count] = h_0[i + on_one - k];
        temp_1[count] = h_1[i + on_one - k];
        temp_2[count] = h_2[i + on_one];
        temp_3[count] = temp_4[count] = 0.0;
        count++;
      }
    }
  }
  for (i = 0; i < on_uc; i++)
  {
    if (onunitcircle[i].imag() < 0)
      k++;
    else if (onunitcircle[i].imag() > 0)
    {
      h_1[i + on_one + on_axis - k] =
          -2 * onunitcircle[i].real() /
          sqrt(onunitcircle[i].real() * onunitcircle[i].real() +
               onunitcircle[i].imag() * onunitcircle[i].imag());
      temp_0[count] = h_0[i + on_one + on_axis - k];
      temp_1[count] = h_1[i + on_one + on_axis - k];
      temp_2[count] = h_2[i + on_one + on_axis];
      temp_3[count] = temp_4[count] = 0.0;
      count++;
    }
  }
  for (i = 0; i < in_uc; i++)
  {
    double r;
    r = sqrt(inside[i].real() * inside[i].real() + inside[i].imag() * inside[i].imag());
    h_2[i + on_one + on_axis + on_uc - k] =
        r * r + 1 / (r * r) +
        4 * inside[i].real() * inside[i].real() /
            (inside[i].real() * inside[i].real() + inside[i].imag() * inside[i].imag());
    h_1[i + on_one + on_axis + on_uc - k] =
        -2 * (r + 1 / r) * inside[i].real() /
        sqrt(inside[i].real() * inside[i].real() + inside[i].imag() * inside[i].imag());
    temp_0[count] = h_0[i + on_one + on_axis + on_uc - k];
    temp_1[count] = h_1[i + on_one + on_axis + on_uc - k];
    temp_2[count] = h_2[i + on_one + on_axis + on_uc - k];
    temp_3[count] = h_1[i + on_one + on_axis + on_uc - k];
    temp_4[count] = 1;
    count++;
  }
  fprintf(file_ptr, "%d\n", count);
  fprintf(file_ptr, "5\n");
  fprintf(file_ptr, "[\n");
  for (i = 0; i < count; i++)
  {
    fprintf(file_ptr, "[%.8f, %.8f, %.8f, %.8f, %.8f],\n", temp_0[i], temp_1[i],
            temp_2[i], temp_3[i], temp_4[i]);
  }
  fprintf(file_ptr, "]\n");
  (void)fclose(file_ptr);
  // save to human readable file
  file_ptr = fopen("FIR_subfilter_response.txt", "w");
  fprintf(file_ptr, " The Subfilters of FIR Filter\n\n");
  fprintf(file_ptr, " The Number of the Subfilters = %d\n", count);
  fprintf(file_ptr, " The Number of impulse response coefficients per subfilter = 5\n\n");
  fprintf(file_ptr, " m hm(0) hm(1) hm(2) hm(3) hm(4)\n");
  fprintf(file_ptr, "====================================================================\n");
  for (i = 0; i < count; i++)
  {
    sig0 = (temp_0[i] >= 0) ? ' ' : '-';
    sig1 = (temp_1[i] >= 0) ? ' ' : '-';
    sig2 = (temp_2[i] >= 0) ? ' ' : '-';
    sig3 = (temp_3[i] >= 0) ? ' ' : '-';
    sig4 = (temp_4[i] >= 0) ? ' ' : '-';
    fprintf(file_ptr, "%3d %c%.8f %c%.8f %c%.8f %c%.8f %c%.8f\n", i + 1, sig0,
            fabs(temp_0[i]), sig1, fabs(temp_1[i]), sig2, fabs(temp_2[i]), sig3,
            fabs(temp_3[i]), sig4, fabs(temp_4[i]));
  }
  (void)fclose(file_ptr);
  delete[] inside;

  delete[] outside;
  delete[] onunitcircle;
  delete[] onaxis;
  delete[] onone;
  delete[] temp_0;
  delete[] temp_1;
  delete[] temp_2;
  delete[] temp_3;
  delete[] temp_4;
  delete[] h_0;
  delete[] h_1;
  delete[] h_2;
  delete[] h_3;
  delete[] h_4;
}
// write the result data file name to monitor screen
void Polynomial ::write_to_screen()
{
  printf(
      "\n\nProgram finished. The results were saved in the follow files:\n\n");
  printf("roots_FIR.dat --------- Roots ang Gain of the FIR filter\n");
  printf("roots_FIR.txt --------- readable text type file\n");
  printf("FIR_subfilter_response.dat --------- subfilters and Gain\n");
  printf("FIR_subfilter_response.txt --------- readable text type file\n\n");
}
// write error message
void Polynomial ::write_error_message(unsigned char error)
{
  printf("Error %d occured!\n", (int)error);
  switch (error)
  {
  case 1:
    printf("Power of polynomial lower null!\n");
    break;
  case 2:
    printf("Polynomial is a null vector!\n");
    break;
  case 3:
    printf("Polynomial is a constant unequal null!\n");
    break;
  }
}
// get the roots we want
unsigned char Polynomial ::get_roots()
{
  const double DBL_EPSILON = 2.2204460492503131E-16;
  std::complex<double> ns; // root determined by Muller's method
  int i;                   // counter
  double newerr;
  unsigned char error; // indicates an error in poly_check
  int red,
      diff; // number of roots at 0
  n -= 1;
  nred = n; // At the beginning: degree defl. polyn. =
  // degree of original polyn.
  maxerr = 0.;
  // check input of the polynomial and make some changes if there are "0" in
  // inputs
  error = poly_check();
  diff = (n - nred); // reduce polynomial, if roots at 0
  p += diff;         // the pointer should change
  n = nred;
  // some errors such like all inputs are Null or "0"
  if (error)
    return error;
  // speical case,polynomial is linear or quadratic,
  // such like ax+b=0 or ax^2 + bx + c=0
  // we can find the result directly and don't need to use Muller & Newton
  // Method
  if (lin_or_quad(p) == 0)
  {
    n += diff; // remember roots at 0
    maxerr = DBL_EPSILON;
    return 0;
  }
  
  std::cout << "INITIAL COEFFS" << std::endl;
  for (int k = 0; k <= nred; k++)
  {
    std::cout << k << " " << p[k] << std::endl;
  }
  monic();
  // Prepare for the input of Muller
  for (i = 0; i <= n; i++)
    pred[i] = p[i];
  
  
  do
  {
    std::cout << " start of loop for nred = " << nred << std::endl;
    for (int k = 0; k <= nred; k++)
    {
      std::cout << k << " " << pred[k] << std::endl;
    }

    // Muller method
    ns = muller();
    // Newton method
    root[nred - 1] = newton(ns, &newerr);
    if (newerr > maxerr)
      maxerr = newerr;
    red = poldef(flag);
    pred += red; // forget lowest coefficients
    nred -= red; // reduce degree of polynomial
    std::cout << " roots of loop for nred = " << nred << std::endl;
    for (int k = 0; k < n; k++)
    {
      std::cout << k << " " << root[k] << std::endl;
    }

  } while (nred > 2);
  // last one or two roots
  (void)lin_or_quad(pred);
  if (nred == 2)
  {
    if (abs(root[1]) <= 1)
    {
      root[1] = newton(root[1], &newerr);
      if (newerr > maxerr)
        maxerr = newerr;
    }
  }
  if (abs(root[0]) <= 1)
    root[0] = newton(root[0], &newerr);
  n += diff; // remember roots at 0
  
  std::cout << "final roots" << std::endl;
  for (int k = 0; k < n; k++)
  {
    std::cout << k << " " << root[k] << std::endl;
  }
  
  if (maxerr < 9e-5)
  {
    printf("\n...\n");
    return 0;
  }
  else
  {
    printf(" Root finding failed, program will exit ...\n");
    exit(0);
  }


  return 0;
}

// monic() computes monic polynomial for original polynomial
void Polynomial ::monic()
{
  double factor; // stores absolute value of the coefficient
  // with highest exponent
  int i;                   // counter variable
  factor = 1. / abs(p[n]); // factor = |1/pn|
  if (factor != 1.)        // get monic pol., when |pn| != 1
    for (i = 0; i <= n; i++)
      p[i] *= factor;
}
// poly_check() check the formal correctness of input
unsigned char Polynomial ::poly_check()
{
  int i = -1, j;
  unsigned char notfound = 1;
  // degree of polynomial less than zero,return error
  if (n < 0)
    return 1;
  // ex. sometimes the degree is 5, but the polynomial is "0,0,3,4,2",
  // so its degree shoule be 3
  for (j = 0; j <= n; j++)
  {
    if (abs(p[j]) != 0.)
      i = j;
  }
  // olynomial is a null
  if (i == -1)
    return 2;
  // polynomials are all "0"
  if (i == 0)
    return 3;
  // get new exponent of polynomial
  n = i;
  i = 0;
  // i --> how many "0" in the input exponent polynomial
  do
  {
    if (abs(p[i]) == 0.)
      i++;
    else
      notfound = 0; // FALSE
  } while (i <= n && notfound);
  if (i == 0)
  { // no '0',original degree=deflated degree
    nred = n;
    return 0;
  }
  else
  { // there are '0', store roots at 0
    for (j = 0; j <= i - 1; j++)
      root[n - j - 1] = std::complex<double>(0., 0.);
    nred = n - i; // reduce degree of deflated polynomial
    return 0;
  }
}
// calculates the roots of a quadratic polynomial ax^2+bx+c=0
void Polynomial ::quadratic(std::complex<double> *p)
{
  std::complex<double> discr, // discriminate
      Z1, Z2,                 // numerators of the quadratic formula
      N;                      // denominator of the quadratic formula
  // discr = p1^2-4*p2*p0
  discr = p[1] * p[1] - 4. * p[2] * p[0];
  // Z1 = -p1+sqrt(discr)
  Z1 = -p[1] + sqrt(discr);
  // Z2 = -p1-sqrt(discr)
  Z2 = -p[1] - sqrt(discr);
  // N = 2*p2
  N = 2. * p[2];
  root[0] = Z1 / N;
  root[1] = Z2 / N;
}

// lin_or_quad() calculates roots of lin. or quadratic equation
unsigned char
Polynomial ::lin_or_quad(std::complex<double> *p)
{
  if (nred == 1)
  { // root = -p0/p1
    root[0] = -p[0] / p[1];
    return 0; // and return no error
  }
  else if (nred == 2)
  { // quadratic polynomial
    quadratic(p);
    return 0; // return no error
  }
  return 1; // nred>2 => no roots were calculated
}
// Horner method to deflate one root
void Polynomial ::hornc(std::complex<double> x0, unsigned char flag)
{
  int i;
  std::complex<double> help1; // help variable
  if ((flag & 1) == 0)        // real coefficients
    for (i = nred - 1; i > 0; i--)
      pred[i].real(pred[i].real() + (x0.real() * pred[i + 1].real()));
  else // complex coefficients
    for (i = nred - 1; i > 0; i--)
    {
      help1 = pred[i + 1] * x0;
      pred[i] = help1 + pred[i];
    }
}
// Horner method to deflate two roots
void Polynomial ::horncd(double a, double b)
{
  int i;
  pred[nred - 1].real(pred[nred - 1].real() + pred[nred].real() * a);
  for (i = nred - 2; i > 1; i--)
  {
    std::cout << i << std::endl;
    pred[i].real(pred[i].real() + (a * pred[i + 1].real() + b * pred[i + 2].real()));
  }
}
// main routine to deflate polynomial
int Polynomial ::poldef(unsigned char flag)
{
  double a, b;
  std::complex<double> x0; // root to be deflated
  x0 = root[nred - 1];
  if (x0.imag() != 0.) // x0 is complex
    flag |= 2;
  if (flag == 2)
  {                    // real coefficients and complex root
    a = 2 * x0.real(); // => deflate x0 and Conjg(x0)
    b = -(x0.real() * x0.real() + x0.imag() * x0.imag());
    root[nred - 2] = conj(x0); // store second root = Conjg(x0)
    horncd(a, b);
    return 2; // two roots deflated
  }
  else
  {
    hornc(x0, flag); // deflate only one root
    return 1;
  }
}
// Newton's method
std::complex<double> Polynomial::newton(std::complex<double> ns, double *dxabs)
{
  const int ITERMAX_1 = 20;
  const double DBL_EPSILON = 2.2204460492503131E-16;
  const double BOUND = sqrt(DBL_EPSILON);
  // if the imaginary part of the root is smaller than BOUND => real root
  const int NOISEMAX = 5;
  const int FACTOR = 5;
  const double FVALUE = 1E36;
  double fabsmin = FVALUE, eps = DBL_EPSILON;
  std::complex<double> x0, // iteration variable for x-value
      xmin,                // best x determined in newton()
      f,                   // P(x0)
      df,                  // P'(x0)
      dx,                  // P(x0)/P'(x0)
      dxh;                 // temperary variable dxh = P(x0)/P'(x0)
  int noise = 0;
  x0 = ns;                           // initial estimation = from Muller method
  xmin = x0;                         // initial estimation for the best x-value
  dx = std::complex<double>(1., 0.); // initial value: P(x0)/P'(x0)=1+j*0
  *dxabs = abs(dx);
  for (iter = 0; iter < ITERMAX_1; iter++)
  {
    f_value1(p, &f, &df, x0); // f=P(x0), df=P'(x0)
    if (abs(f) < fabsmin)
    {
      xmin = x0;
      fabsmin = abs(f);
      noise = 0;
    }
    if (abs(df) != 0.)
    { // calculate new dx
      dxh = f / df;
      if (abs(dxh) < *dxabs * FACTOR)
      {
        dx = dxh;
        *dxabs = abs(dx);
      }
    }
    if (abs(xmin) != 0.)
    {
      if (*dxabs / abs(xmin) < eps || noise == NOISEMAX)
      {
        if (fabs(xmin.imag()) < BOUND && flag == 0)
        {
          xmin.imag(0.); // if imag. part<BOUND, let's it=0
        }
        *dxabs = *dxabs / abs(xmin);
        return xmin; // return best approximation
      }
    }
    // x0 = x0 - P(x0)/P'(x0)
    x0 -= dx;
    noise++;
  }
  if (fabs(xmin.imag()) < BOUND && flag == 0)
    xmin.imag(0.);
  // if imag. part<BOUND , let's it=0
  if (abs(xmin) != 0.)
    *dxabs = *dxabs / abs(xmin);
  return xmin; // return best xmin until now
}
void Polynomial ::f_value1(std::complex<double> *p, std::complex<double> *f, std::complex<double> *df, std::complex<double> x0)
{
  int i;                      // counter
  std::complex<double> help1; // temperary variable
  *f = p[n];
  // COMPLEXM(*df, 0., 0.);
  df->real(0.);
  df->imag(0.);
  for (i = n - 1; i >= 0; i--)
  {
    *df = (*df) * x0 + (*f);
    *f = (*f) * x0 + p[i];
  }
}
void Polynomial ::f_value2(std::complex<double> *f, std::complex<double> *df, std::complex<double> x0)
{
  int i;                      // counter
  std::complex<double> help1; // temperary variable
  *f = psave[nred];
  // COMPLEXM(*df, 0., 0.);
  df->real(0.);
  df->imag(0.);
  for (i = nred - 1; i >= 0; i--)
  {
    *df = (*df) * x0 + (*f);
    *f = (*f) * x0 + psave[i];
  }
  
}
// Muller's method
std::complex<double> Polynomial ::muller()
{
  const int ITERMAX = 150;    // max. number of iteration steps
  const double FVALUE = 1e36; // initialisation of |P(x)|^2
  const double DBL_EPSILON = 2.2204460492503131E-16;
  const double NOISESTART = DBL_EPSILON * 1e2;
  const int NOISEMAX = 5;
  double h2abs,         // h2abs=|h2| h2absnew=distance between old and new x2
      f1absq,           // f1absq=|f1|^2 used for check
      f2absq = FVALUE,  // f2absq=|f2|^2 used for check
      f2absqb = FVALUE, // f2absqb=|P(xb)|^2 used for check
      epsilon;
  int seconditer = 0, // second iteration, when root is too bad
      noise = 0,      // noise counter
      rootd = 0;
  std::complex<double> xb;   // best x-value
  initialize(&xb, &epsilon); // initialize x0,x1,x2,h1,h2,q2,*xb
  // use Horner's Method, get f0=P(x0), f1=P(x1), f2=P(x2)
  f_value(nred, &f0, x0);
  f_value(nred, &f1, x1);
  f_value(nred, &f2, x2);
  do
  {
    do
    {
      // get q2 ( q2=2C/B(+-)sqr(B^2-4AC) )
      root_of_parabola();
      // store values for the next iteration
      x0 = x1;
      x1 = x2;
      h2abs = abs(h2); // |x2-x1|
                       // get the result from Muller's method: x2=x2-(x2-x1) *
                       // 2C/B(+-)sqr(B^2-4AC)
      iteration_equation(&h2abs);
      // store P(x) values for the next iteration
      f0 = f1;
      f1 = f2;
      f1absq = f2absq;
      compute_function(f1absq, &f2absq, epsilon);
      // check if the new x2 is best enough , these two checks are necessary
      check_x_value(&xb, &f2absqb, &rootd, f1absq, f2absq, epsilon, &noise);
      if (fabs((abs(xb) - abs(x2)) / abs(xb)) < NOISESTART)
        noise++;
    } while ((iter <= ITERMAX) && (!rootd) && (noise <= NOISEMAX));
    seconditer++;
    root_check(f2absqb, &seconditer, &rootd, &noise, xb);
  } while (seconditer == 2);
  return xb; // return best x value
}
// initializing routine
void Polynomial ::initialize(std::complex<double> *xb, double *epsilon)
{
  const double DBL_EPSILON = 2.2204460492503131E-16;
  x0 = std::complex<double>(0., 1.);                     // x0 = 0 + j*1
  x1 = std::complex<double>(0., -1.);                    // x1 = 0 - j*0
  x2 = std::complex<double>(1. / sqrt(2), 1. / sqrt(2)); // x2 = (1 + j*1)/sqrt(2)
  h1 = x1 - x0;
  h2 = x2 - x1; // h2 = x2 - x1
  q2 = h2 / h1; // q2 = h2/h1
  *xb = x2;     // best initial x-value = zero
  *epsilon = 5 * DBL_EPSILON;
  iter = 0; // reset iteration counter
}
// root of Muller's parabola------q2
void Polynomial ::root_of_parabola(void)
{
  std::complex<double> A2, B2, C2, discr, N1, N2;
  const double DBL_EPSILON = 2.2204460492503131E-16;
  // A2 = q2(f2 - (1+q2)f1 + f0q2)
  // B2 = q2[q2(f0-f1) + 2(f2-f1)] + (f2-f1)
  // C2 = (1+q2)f[2]
  A2 = q2 * (f2 - (1. + q2) * f1 + f0 * q2);
  B2 = q2 * (q2 * (f0 - f1) + 2. * (f2 - f1)) + (f2 - f1);
  C2 = (1. + q2) * f2;
  // discr = B2^2 - 4A2C2
  discr = B2 * B2 - 4. * A2 * C2;
  // denominators of q2
  N1 = B2 - sqrt(discr);
  N2 = B2 + sqrt(discr);
  // choose denominater with largest modulus
  if (abs(N1) > abs(N2) && abs(N1) > DBL_EPSILON)
    q2 = (-2.) * C2 / N1;
  else if (abs(N2) > DBL_EPSILON)
    q2 = (-2.) * C2 / N2;
  else
    q2 = std::complex<double>(cos(iter), sin(iter));
}
// main iteration equation: x2 = h2*q2 + x2
void Polynomial ::iteration_equation(double *h2abs)
{
  double h2absnew, // Absolute value of the new h2
      help;        // help variable
  const double MAXDIST = 1e3;
  h2 *= q2;
  h2absnew = abs(h2); // distance between old and new x2
  if (h2absnew > (*h2abs * MAXDIST))
  { // maximum relative change
    help = MAXDIST / h2absnew;
    h2 *= help;
    q2 *= help;
  }
  *h2abs = h2absnew; // actualize old distance for next iteration
  x2 += h2;
}
// use Horner's method to get P(x)
void Polynomial ::f_value(int n, std::complex<double> *f, std::complex<double> x0)
{
  int i;
  std::complex<double> help1;
  *f = pred[n];
  // compute P(x0)
  for (i = n - 1; i >= 0; i--)
  {
    // use Horner's method
    help1 = *f * x0; // *f = p[i] + *f * x0
    *f = help1 + pred[i];
  }
}
// check of too big function values
void Polynomial ::too_big_functionvalues(double *f2absq)
{
  const double DBL_MAX = 1.7976931348623157E+308;
  const double BOUND4 = sqrt(DBL_MAX) / 1e4;
  if ((fabs(f2.real()) + fabs(f2.imag())) > BOUND4) // limit |f2|^2
    *f2absq = fabs(f2.real()) + fabs(f2.imag());
  else
    *f2absq = (f2.real()) * (f2.real()) + (f2.imag()) * (f2.imag());
}
void Polynomial::suppress_overflow()
{
  int kiter;          // internal iteration counter
  unsigned char loop; // loop = FALSE => terminate loop
  double help;        // help variable
  const double KITERMAX = 1e3;
  const double DBL_MAX = 1.7976931348623157E+308;
  const double BOUND4 = sqrt(DBL_MAX) / 1e4;
  const double BOUND6 = log10(BOUND4) - 4;
  kiter = 0; // reset iteration counter
  do
  {
    loop = 0;       // initial estimation: no overflow
    help = abs(x2); // help = |x2|
    if (help > 1. && fabs(nred * log10(help)) > BOUND6)
    {
      kiter++; // if |x2|>1 and |x2|^nred>10^BOUND6
      if (kiter < KITERMAX)
      {               // then halve the distance between
        h2 = .5 * h2; // new and old x2
        q2 = .5 * q2;
        x2 = x2 - h2;
        loop = 1;
      }
      else
        kiter = 0;
    }
  } while (loop);
}
// Muller's modification to improve convergence
void Polynomial::convergence_check(int *overflow, double f1absq, double f2absq,
                                   double epsilon)
{
  const int CONVERGENCE = 100;
  const int ITERMAX = 150;
  if ((f2absq > (CONVERGENCE * f1absq)) && (abs(q2) > epsilon) &&
      (iter < ITERMAX))
  {
    q2 *= .5; // in case of overflow:
    h2 *= .5; // halve q2 and h2; compute new x2
    x2 -= h2;
    *overflow = 1;
  }
}
// compute P(x2) and make some checks
void Polynomial ::compute_function(double f1absq, double *f2absq,
                                   double epsilon)
{
  int overflow; // overflow = TRUE => overflow occures
  // overflow = FALSE => no overflow occures
  do
  {
    overflow = 0; // initial estimation: no overflow
    // suppress overflow
    suppress_overflow();
    // calculate new value => result in f2
    f_value(nred, &f2, x2);
    // check of too big function values
    too_big_functionvalues(f2absq);
    // increase iterationcounter
    iter++;
    // Muller's modification to improve convergence
    convergence_check(&overflow, f1absq, *f2absq, epsilon);
  } while (overflow);
}
// check if the new x2 the best approximation
void Polynomial ::check_x_value(std::complex<double> *xb, double *f2absqb, int *rootd,
                                double f1absq, double f2absq, double epsilon,
                                int *noise)
{
  const double BOUND1 = 1.01;
  const double BOUND2 = 0.99;
  const double BOUND3 = 0.01;
  if ((f2absq <= (BOUND1 * f1absq)) && (f2absq >= (BOUND2 * f1absq)))
  {
    // function-value changes slowly
    if (abs(h2) < BOUND3)
    {          // if |h[2]| is small enough =>
      q2 *= 2; // double q2 and h[2]
      h2 *= 2;
    }
    else
    { // otherwise: |q2| = 1 and
      // h[2] = h[2]*q2
      q2 = std::complex<double>(cos(iter), sin(iter));
      h2 = h2 * q2;
    }
  }
  else if (f2absq < *f2absqb)
  {
    *f2absqb = f2absq; // the new function value is the
    *xb = x2;          // best approximation
    *noise = 0;        // reset noise counter
    if ((sqrt(f2absq) < epsilon) && (abs((x2 - x1) / x2)) < epsilon)
      *rootd = 1;
  }
}
void Polynomial ::root_check(double f2absqb, int *seconditer, int *rootd,
                             int *noise, std::complex<double> xb)
{
  std::complex<double> df; // df=P'(x0)
  const double BOUND7 = 1e-5;
  if ((*seconditer == 1) && (f2absqb > 0))
  {
    f_value2(&f2, &df, xb); // f2=P(x0), df=P'(x0)
    if (abs(f2) / (abs(df) * abs(xb)) > BOUND7)
    {
      // start second iteration with new initial estimations
      x0 = std::complex<double>(-1. / sqrt(2), 1. / sqrt(2));
      x1 = std::complex<double>(1. / sqrt(2), -1. / sqrt(2));
      x2 = std::complex<double>(-1. / sqrt(2), -1. / sqrt(2));
      f_value(nred, &f0, x0);
      f_value(nred, &f1, x1);
      f_value(nred, &f2, x2);
      iter = 0;        // reset iteration counter
      (*seconditer)++; // increase seconditer
      *rootd = 0;      // no root determined
      *noise = 0;      // reset noise counter
    }
  }
}
void Polynomial ::check_value()
{
  int k, m, i;
  double omega, W_Re, W_Im, T_Re, T_Im, H_mag, upper, lower, G;
  double s_H_mag;
  double s_T_Re, s_T_Im, s_H_Re, s_H_Im, s_temp;
  std::complex<double> H, temp;
  double *test2;
  double *h;
  test2 = new double[NPLOT];
  FILE *file_ptr;
  H.real(1);
  H.imag(0);
  upper = lower = 0.0;
  h = new double[N + 1];
  for (i = N - 1; i >= 0; i--)
    h[i] = matlab[i];
  for (k = 0; k < NPLOT; k++)
  {
    omega = PI * ((double)k) / ((double)NPLOT);
    W_Re = cos(omega);
    W_Im = -sin(omega);
    T_Re = cos(2 * omega);
    T_Im = -sin(2 * omega);
    s_T_Re = 1.0;
    s_T_Im = 0.0;
    s_H_Re = h[0];
    s_H_Im = 0.0;
    for (m = 0; m < N - 1; m++)
    {
      if (root[m].imag() == 0.0)
      {
        temp.real(1 - root[m].real() * W_Re);
        temp.imag(root[m].real() * W_Im);
      }
      else
      {
        temp.real(1 - (2 * root[m].real() * W_Re) +
                  (root[m].real() * root[m].real() + (root[m].imag() * root[m].imag())) * T_Re);
        temp.imag((root[m].real() * root[m].real() + (root[m].imag() * root[m].imag())) * T_Im -
                  2 * root[m].real() * W_Im);
        m++;
      }
      H *= temp;
    }
    for (i = 1; i < N; i++)
    {
      s_temp = s_T_Re * W_Re - s_T_Im * W_Im;
      s_T_Im = s_T_Re * W_Im + s_T_Im * W_Re;
      s_T_Re = s_temp;
      if (i < distinct)
      {
        s_temp = h[i];
      }
      else
      {
        s_temp = h[N - 1 - i];
        if (indicator == 1)
          s_temp = -s_temp;
      }
      s_H_Re += s_temp * s_T_Re;
      s_H_Im += s_temp * s_T_Im;
    }
    s_H_mag = sqrt(s_H_Re * s_H_Re + s_H_Im * s_H_Im);
    H_mag = sqrt(H.real() * H.real() + H.imag() * H.imag());
    H.real(1);
    H.imag(0);
    upper += s_H_mag * H_mag;
    lower += H_mag * H_mag;
    test2[k] = s_H_mag;
  }
  G = upper / lower;
  file_ptr = fopen("roots_FIR.dat", "a");
  fprintf(file_ptr, "%.10f\n", G);
  fclose(file_ptr);
  file_ptr = fopen("FIR_subfilter_response.dat", "a");
  fprintf(file_ptr, "%.10f\n", G);
  fclose(file_ptr);
  file_ptr = fopen("roots_FIR.txt", "a");
  fprintf(file_ptr, "Gain = %.10f\n", G);
  fclose(file_ptr);
  file_ptr = fopen("FIR_subfilter_response.txt", "a");
  fprintf(file_ptr, "Gain = %.10f\n", G);
  fclose(file_ptr);
  delete[] matlab;
  delete[] h;
  delete[] test2;
}
Polynomial ::~Polynomial()
{
  delete[] psave1;
  delete[] root;
  delete[] psave;
}
int main(void)
{
  Polynomial Poly;
  unsigned char error;
  Poly.read_from_file();
  error = Poly.get_roots();
  if (!error)
  {
    Poly.write_to_file();
    Poly.write_to_screen();
    Poly.check_value();
  }
  else
  {
    Poly.write_error_message(error);
    return 1;
  }
  return 0;
}