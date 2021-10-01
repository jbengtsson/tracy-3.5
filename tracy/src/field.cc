 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */


#if NO_TPSA == 1

void tps::print(const string &str)
{
  std::cout << std::scientific << std::setprecision(6) << str << "cst:\n"
	    << std::setw(14) << ltps[0] << "\nlinear:\n";
  for (int i = 1; i <= nv_tps; i++)
    std::cout << std::scientific << std::setprecision(6)
	      << std::setw(14) << ltps[i];
  std::cout << "\n";
}

template<>
void ss_vect<tps>::print(const string &str)
{
  std::cout << str;
  for (int i = 1; i <= nv_tps; i++) {
    for (int j = 1; j <= nv_tps; j++)
      std::cout << std::scientific << std::setprecision(6)
		<< std::setw(14) << get_m_ij(*this, i, j);
    std::cout << "\n";
  }
}

#else

void tps::print(const string &str) { std::cout << str << *this; }

template<>
void ss_vect<tps>::print(const string &str) { std::cout << str << *this; }

#endif

template<>
void ss_vect<double>::print(const string &str)
{
  std::cout << std::scientific << std::setprecision(6) << str << std::setw(14)
	    << *this << "\n"; 
}


template<typename T>
ss_vect<T>& ss_vect<T>::operator=(const ss_vect<T> &x)
{
  for (int i = 0; i < ps_dim; i++)
    (*this)[i] = x[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator+=(const ss_vect<T> &a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] += a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator-=(const ss_vect<T> &a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] -= a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator*=(const double a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] *= a;
  return *this;
}

template<>
ss_vect<tps>& ss_vect<tps>::operator*=(const tps &a)
{
  for (int i = 0; i < ps_dim; i++)
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
  for (int i = 0; i < ps_dim; i++)
    ss[i] = 0e0;
  return *this;
}


template<>
ss_vect<double> ss_vect<double>::identity(void)
{
  printf("\nidentity: not implemented for ss_vect<double>\n");
  exit(1);
}

template<>
ss_vect<tps> ss_vect<tps>::identity(void)
{
  for (int i = 0; i < ps_dim; i++)
   ss[i] = tps(0e0, i+1);
  return *this;
}


template<typename CharT, class Traits>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, ss_vect<tps> &a)
{
  tps b;

  for (int i = 0; i < ps_dim; i++)
    is >> a[i];
  return is;
}

// instantiate
template std::basic_istream<char, std::char_traits<char> >&
operator>>(std::basic_istream<char, std::char_traits<char> > &, ss_vect<tps> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<double> &a)
{
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (int i = 0; i < ps_dim; i++)
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
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (int i = 0; i < ps_dim; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<tps> &);


