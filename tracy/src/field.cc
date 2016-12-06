 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

template<>
ss_vect<double> ss_vect<tps>::cst(void) const
{
  int              i;
  ss_vect<double>  x;

  for (i = 0; i < ss_dim; i++)
    x[i] = (*this)[i].cst();
  return x;
}


template<typename T>
ss_vect<T>& ss_vect<T>::operator=(const ss_vect<T> &x)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    (*this)[i] = x[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator+=(const ss_vect<T> &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] += a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator-=(const ss_vect<T> &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] -= a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator*=(const double a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] *= a;
  return *this;
}

template<>
ss_vect<tps>& ss_vect<tps>::operator*=(const tps &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] *= a;
  return *this;
}


template<typename T>
ss_vect<T> operator+(const ss_vect<T> &x) { return ss_vect<T>(x); }

template<typename T>
ss_vect<T> operator-(const ss_vect<T> &x) { return ss_vect<T>(x) *= -1; }

// instantiate
template ss_vect<double> operator-(const ss_vect<double> &);
template ss_vect<tps> operator-(const ss_vect<tps> &);

ss_vect<double> operator+(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<double> operator-(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<double> operator*(const ss_vect<double> &a, const double b)
{ return ss_vect<double>(a) *= b; }

ss_vect<double> operator*(const double a, const ss_vect<double> &b)
{ return ss_vect<double>(b) *= a; }

ss_vect<tps> operator*(const ss_vect<tps> &a, const double b)
{ return ss_vect<tps>(a) *= b; }

ss_vect<tps> operator*(const double a, const ss_vect<tps> &b)
{ return ss_vect<tps>(b) *= a; }


template<typename T>
ss_vect<T> ss_vect<T>::zero(void)
{
  int         i;

  for (i = 0; i < ss_dim; i++)
    ss[i] = 0.0;
  return *this;
}

template<>
ss_vect<tps> ss_vect<tps>::identity(void)
{
  int           i;

  for (i = 0; i < ss_dim; i++)
   ss[i] = tps(0.0, i+1);
  return *this;
}


template<typename CharT, class Traits>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, ss_vect<tps> &a)
{
  int  i;
  tps  b;

  for (i = 0; i < 6; i++) {
    is >> a[i];
  }
  return is;
}

// instantiate
template std::basic_istream<char, std::char_traits<char> >&
operator>>(std::basic_istream<char, std::char_traits<char> > &, ss_vect<tps> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<double> &a)
{
  int                                 i;
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (i = 0; i < 6; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
//   s << endl;
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<double> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<tps> &a)
{
  int                                 i;
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (i = 0; i < 6; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<tps> &);


